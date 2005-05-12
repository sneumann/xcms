void NetCDFStrError(const int *ncerr, const int *len, char *errortext[]);

void NetCDFOpen(const char *fileName[], int *ncid, int *status);

void NetCDFClose(const int *ncid, int *status);

void NetCDFVarID(const int *ncid, const char *varName[], int *varid, int *status);

void NetCDFVarLen(const int *ncid, const int *varid, int *len, int *status);

void NetCDFVarDouble(const int *ncid, const int *varid, double *data, int *status);

void NetCDFVarInt(const int *ncid, const int *varid, int *data, int *status);

void NetCDFMSPoints(const int *ncid, const int *scanNumber, 
                    const int *scanIndex, const int *pointNumber, 
                    double *massValues, double *intensityValues, int *status);
