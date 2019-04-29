REM files which must come first
REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 addHHtoTjm.f95
REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 calcFjm.f95
REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 CalcRjm.f95
REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 calctae.f95
REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 calctjm.f95

REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 EvaluateCombination.f95
REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 fnames.f95
REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 initcomb.f95
REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 readpop.f95
REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 ReadTableInfo.f95

REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 readujm.f95
REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 removeHHfromTjm.f95
REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 replace_RSSZm.f95
REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 replace_TAE.f95
REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 Write_Full_weights.f95

REM gfortran -c -Wall -Wextra -Warray-temporaries -Wconversion -fbounds-check -O1 CO.f95

REM Shared objects
REM gfortran -g -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal -o "..\CO Dummy\CO\CO_BU.exe"^
REM    CO.o addHHtoTjm.o calcFjm.o CalcRjm.o calctae.o calctjm.o EvaluateCombination.o^
REM	fnames.o initcomb.o readpop.o ReadTableInfo.o readujm.o removeHHfromTjm.o replace_RSSZm.o replace_TAE.o Write_Full_weights.o^

REM -Warray-temporaries 

gfortran -o "..\CO Dummy\CO\CO_BU.exe" CO.f95 ^
	addHHtoTjm.f95 calcFjm.f95 CalcRjm.f95 calctae.f95 calctjm.f95 EvaluateCombination.f95^
	fnames.f95 initcomb.f95 readpop.f95 ReadTableInfo.f95 readujm.f95 removeHHfromTjm.f95^
	replace_RSSZm.f95 replace_TAE.f95 Write_Full_weights.f95^
	-g -static -Wall -Wextra -Wconversion -fimplicit-none -fbacktrace^
	-ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,denormal,zero,overflow,underflow -finit-real=nan
					

