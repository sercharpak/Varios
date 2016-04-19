all: string.x string_100.000000.dat cuerda_100.0.pdf

cuerda_100.0.pdf:string_100.000000.dat dasWave.py
	python dasWave.py string_100.000000.dat
string_100.000000.dat :string.x
	./string.x 100
string.x: dasPunto1.c
	cc dasPunto1.c -lm -o string.x 
clean:
	rm *.dat *.pdf
