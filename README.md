YOU MUST COMPILE CMS LIKE:

```
git clone https://github.com/msoos/cryptominisat.git
cd cryptominisat
git checkout scalmc-simpsat2
rm -rf build
mkdir build
cd build
cmake -DUSE_GAUSS=ON ..
make -j4
sudo make install
```

If you don't it will NOT work and will NOT compile. If you made a mistake please re-do EVERYTHING starting at "rm -rf build". Yes, you MUST delete the "build" directory. Otherwise you WILL NOT be able to compile.
