@ECHO OFF
SETLOCAL EnableDelayedExpansion

:: Initialize to directory of data
ECHO Setting directory...
set DIR="C:/Users/jvons/Documents/NCF/Thesis/Datasets"
cd /d %DIR%
ECHO Directory set to %DIR%.

:: Grab data file(s)
set datalist=HPNE PANC-1 Capan-1
(for %%a in (%datalist%) do (
	set datadir=%DIR% %%a
	ECHO %datadir%))

PAUSE
