#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include "netCDF.h"

void NetCDFStrError(const int *ncerr, const int *len, char *errortext[]) {

    strncpy(errortext[0], nc_strerror(*ncerr), *len);
}

void NetCDFOpen(const char *fileName[], int *ncid, int *status) {
    
    *status = nc_open(fileName[0], NC_NOWRITE, ncid);
}

void NetCDFClose(const int *ncid, int *status) {
    
    *status = nc_close(*ncid);
}

void NetCDFVarID(const int *ncid, const char *varName[], int *varid, int *status) {
    
    *status = nc_inq_varid(*ncid, varName[0], varid);
}

void NetCDFVarLen(const int *ncid, const int *varid, int *len, int *status) {
    
    int    ndims, dimids[NC_MAX_VAR_DIMS], i;
    size_t dimLen;
    
    if ((*status = nc_inq_varndims(*ncid, *varid, &ndims)))
        return;
    
    if ((*status = nc_inq_vardimid(*ncid, *varid, dimids)))
        return;
    
    *len = 1;
    for (i = 0; i < ndims; i++) {
        if ((*status = nc_inq_dimlen(*ncid, dimids[i], &dimLen)))
            return;
        *len *= dimLen;
    }
}

void NetCDFVarDouble(const int *ncid, const int *varid, double *data, int *status) {
    
    int    varLen, i;
    double scaleFactor, addOffset;
    size_t attLen;
    
    NetCDFVarLen(ncid, varid, &varLen, status);
    if (*status)
        return;
    
    if ((*status = nc_get_var_double(*ncid, *varid, data)))
        return;
    
    if (!nc_inq_att(*ncid, *varid, "scale_factor", NULL, &attLen))
        if (attLen == 1 && !nc_get_att_double(*ncid, *varid, "scale_factor", &scaleFactor) && scaleFactor != 1)
            for (i = 0; i < varLen; i++)
                data[i] *= scaleFactor;
    
    if (!nc_inq_att(*ncid, *varid, "add_offset", NULL, &attLen))
        if (attLen == 1 && !nc_get_att_double(*ncid, *varid, "add_offset", &addOffset) && addOffset != 0)
            for (i = 0; i < varLen; i++)
                data[i] += addOffset;
}

void NetCDFVarInt(const int *ncid, const int *varid, int *data, int *status) {
    
    int    varLen;
    
    NetCDFVarLen(ncid, varid, &varLen, status);
    if (*status)
        return;
    
    *status = nc_get_var_int(*ncid, *varid, data);
}

void NetCDFMSPoints(const int *ncid, const int *scanNumber, 
                    const int *scanIndex, const int *pointNumber, 
                    double *massValues, double *intensityValues, int *status) {
    
    int    varid, i, j, scanLen;
    double tmpMass, tmpIntensity;
    
    *status = nc_inq_varid(*ncid, "mass_values", &varid);
    if (*status)
        return;
    
    NetCDFVarDouble(ncid, &varid, massValues, status);
    if (*status)
        return;
    
    *status = nc_inq_varid(*ncid, "intensity_values", &varid);
    if (*status)
        return;
    
    NetCDFVarDouble(ncid, &varid, intensityValues, status);
    if (*status)
        return;
    
    for (i = 0; i < *scanNumber-1; i++)
        if (scanIndex[i+1] - scanIndex[i] > 1) 
            if (massValues[scanIndex[i]] < massValues[scanIndex[i]+1])
                return;
    
    for (i = 0; i < *scanNumber; i++) {
        scanLen = (i < *scanNumber - 1) ? scanIndex[i+1] - scanIndex[i] :
                                          *pointNumber - scanIndex[i];
        for (j = 0; j < scanLen/2; j++) {
            tmpMass = massValues[scanIndex[i]+j];
            tmpIntensity = intensityValues[scanIndex[i]+j];
            massValues[scanIndex[i]+j] = massValues[scanIndex[i]+scanLen-1-j];
            intensityValues[scanIndex[i]+j] = intensityValues[scanIndex[i]+scanLen-1-j];
            massValues[scanIndex[i]+scanLen-1-j] = tmpMass;
            intensityValues[scanIndex[i]+scanLen-1-j] = tmpIntensity;
        }
    }
}
