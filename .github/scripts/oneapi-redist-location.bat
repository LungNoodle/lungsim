@echo off
for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
echo "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\windows\redist\intel64_win\compiler"

