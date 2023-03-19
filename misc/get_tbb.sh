wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz && rm 2019_U9.tar.gz
mkdir -p $(dirname -- "$0";)"/../lib"
mv -vn oneTBB-2019_U9 $(dirname -- "$0";)"/../lib/tbb"
make -C "$(dirname -- "$0";)/../lib/tbb"
