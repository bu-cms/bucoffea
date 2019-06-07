XRDURL="http://xrootd.org/download/v4.9.1/xrootd-4.9.1.tar.gz"
# wget $XRDURL
# tar xf $(basename $XRDURL)

XRDBASE=$(readlink -e ./$(echo $(basename $XRDURL) | sed "s|.tar.gz||"))
BUILDDIR=$XRDBASE/build
mkdir -p $BUILDDIR
cd $BUILDDIR




mkdir -p external
pushd external
if [ ! -d xrootd ]; then
    git clone git@github.com:xrootd/xrootd.git
fi

pushd xrootd
mkdir build
pushd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j4
make install -j4

pip install -e bindings/python



