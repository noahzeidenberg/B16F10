@echo off
REM Batch script to run pathway enrichment analysis for multiple GSE IDs

REM Check if GSE ID list file is provided
if "%~1"=="" (
    echo Error: No GSE ID list file provided
    echo Usage: run_pathway_enrichment_batch.bat GSE_ID_LIST_FILE
    exit /b 1
)

REM Get GSE ID list file from command line argument
set GSE_ID_LIST_FILE=%~1

REM Check if the file exists
if not exist "%GSE_ID_LIST_FILE%" (
    echo Error: GSE ID list file %GSE_ID_LIST_FILE% does not exist
    exit /b 1
)

REM Create logs directory if it doesn't exist
if not exist "logs" mkdir logs

REM Read GSE IDs from the file and submit jobs
for /f "tokens=*" %%a in ('type "%GSE_ID_LIST_FILE%"') do (
    echo Submitting job for GSE ID: %%a
    sbatch scripts/run_pathway_enrichment.slurm %%a
    REM Wait a bit to avoid overwhelming the scheduler
    timeout /t 5 /nobreak > nul
)

echo All jobs submitted successfully
exit /b 0 