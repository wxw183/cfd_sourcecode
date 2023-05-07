gcc -o adi.exe adi.c -fopenmp
mkdir 01
copy /y adi.exe 01
copy /y plot3d.py 01
mkdir 02
copy /y adi.exe 02
copy /y plot3d.py 02
mkdir 03
copy /y adi.exe 03
copy /y plot3d.py 03
mkdir 04
copy /y adi.exe 04
copy /y plot3d.py 04
mkdir 05
copy /y adi.exe 05
copy /y plot3d.py 05
mkdir 06
copy /y adi.exe 06
copy /y plot3d.py 06
mkdir 07
copy /y adi.exe 07
copy /y plot3d.py 07
cd 01
adi 0.1 0.02 10 8
python plot3d.py
cd ..
cd 02
adi 0.05 0.02 10 8
python plot3d.py
cd ..
cd 03
adi 0.025 0.02 10 8
python plot3d.py
cd ..
cd 04
adi 0.05 0.005 10 8
python plot3d.py
cd ..
cd 05
adi 0.05 0.01 10 8
python plot3d.py
cd ..
cd 06
adi 0.05 0.02 10 8
python plot3d.py
cd ..
cd 07
adi 0.05 0.04 10 8
python plot3d.py
cd ..
python adi.py
