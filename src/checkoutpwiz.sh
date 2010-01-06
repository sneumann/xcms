
#rm -rf pwiz
#mkdir pwiz
cd pwiz

PWIZREPO=https://proteowizard.svn.sourceforge.net/svnroot/proteowizard/trunk/pwiz/pwiz

#for DIR in data/msdata data/common utility/misc/ utility/minimxml/
#for DIR in  ; do 
#    svn co $PWIZREPO/$DIR $DIR
#done

#------------------

cd ..
#mkdir boost
cd boost


BOOSTVER=Boost_1_39_0
BOOSTREPO=http://svn.boost.org/svn/boost/tags/release/$BOOSTVER/boost

#svn co --non-recursive $BOOSTREPO .

#for DIR in smart_ptr  config config mpl detail iostreams exception  \
# type_traits preprocessor format algorithm logic optional range \
# iterator function utility concept bind regex filesystem system thread \
# lambda  tuple multi_index serialization archive functional integer ; do 
#for DIR in   ; do 
#    svn co $BOOSTREPO/$DIR $DIR
#done

BOOSTREPO=http://svn.boost.org/svn/boost/tags/release/$BOOSTVER/libs
#for DIR in iostreams/src thread/src/pthread/ filesystem/src/ regex/src ; do 
for DIR in system/src ; do 
    svn co $BOOSTREPO/$DIR $DIR
done



