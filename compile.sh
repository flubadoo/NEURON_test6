gcc -v -I/usr/include/python2.7 -L./matrix/ -lpython2.7 -shared -o rxd.so -fPIC rxd.c matrix/meschach.a -lm
