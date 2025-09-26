@echo off

echo. > output.log
set CCES_DIR=C:\analog\cces\3.0.2
set logfile=output.log
set searchstring=Test pass

start "OPENOCD" "%CCES_DIR%\ARM\openocd\bin\openocd.exe" -f interface\ice1000.cfg -f board\adspsc835w_ev_som.cfg

start /wait "XT-GDB" "%CCES_DIR%\Xtensa\XtensaTools\bin\xt-gdb.exe" --batch --command="cmd_audiomark/gdb_commands_audiomark.gdb" -l 50

findstr /c:"%searchstring%" "%logfile%" >nul

taskkill -f /im openocd.exe

start "OPENOCD" "%CCES_DIR%\ARM\openocd\bin\openocd.exe" -f interface\ice1000.cfg -f board\adspsc835w_ev_som.cfg

start /wait "XT-GDB" "%CCES_DIR%\Xtensa\XtensaTools\bin\xt-gdb.exe" --batch --command="cmd_audiomark/gdb_commands_aec.gdb" -l 50

findstr /c:"%searchstring%" "%logfile%" >nul

taskkill -f /im openocd.exe

start "OPENOCD" "%CCES_DIR%\ARM\openocd\bin\openocd.exe" -f interface\ice1000.cfg -f board\adspsc835w_ev_som.cfg

start /wait "XT-GDB" "%CCES_DIR%\Xtensa\XtensaTools\bin\xt-gdb.exe" --batch --command="cmd_audiomark/gdb_commands_abf.gdb" -l 50

findstr /c:"%searchstring%" "%logfile%" >nul

taskkill -f /im openocd.exe

start "OPENOCD" "%CCES_DIR%\ARM\openocd\bin\openocd.exe" -f interface\ice1000.cfg -f board\adspsc835w_ev_som.cfg

start /wait "XT-GDB" "%CCES_DIR%\Xtensa\XtensaTools\bin\xt-gdb.exe" --batch --command="cmd_audiomark/gdb_commands_anr.gdb" -l 50

findstr /c:"%searchstring%" "%logfile%" >nul

taskkill -f /im openocd.exe

start "OPENOCD" "%CCES_DIR%\ARM\openocd\bin\openocd.exe" -f interface\ice1000.cfg -f board\adspsc835w_ev_som.cfg

start /wait "XT-GDB" "%CCES_DIR%\Xtensa\XtensaTools\bin\xt-gdb.exe" --batch --command="cmd_audiomark/gdb_commands_mfcc.gdb" -l 50

findstr /c:"%searchstring%" "%logfile%" >nul

taskkill -f /im openocd.exe

start "OPENOCD" "%CCES_DIR%\ARM\openocd\bin\openocd.exe" -f interface\ice1000.cfg -f board\adspsc835w_ev_som.cfg

start /wait "XT-GDB" "%CCES_DIR%\Xtensa\XtensaTools\bin\xt-gdb.exe" --batch --command="cmd_audiomark/gdb_commands_kws.gdb" -l 50

findstr /c:"%searchstring%" "%logfile%" >nul

taskkill -f /im openocd.exe


findstr /v /c:"logging" ^
              /c:"logging enabled:  on: Logging is enabled." ^
              /c:"warning:" ^
              /c:"Loading" ^
              /c:"Program" ^
              /c:"determining" ^
 		/c:"_exit ()" ^
		/c:"83" ^
		/c:"[" ^
		/c:"0x20211" ^
              /c:"Loading section" output.log > final_output.log 2>nul

echo complete. Output result saved to final_output.log.
