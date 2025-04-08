@echo off
setlocal enabledelayedexpansion

echo Starting design matrix generation...

:: Set absolute paths - EDIT THESE PATHS
set LM_STUDIO_PATH="C:\Users\nbfz0\.lmstudio\bin\lms.exe"
set MODEL_PATH="gemma-2-9b-it-GGUF/gemma-2-9b-it-Q4_K_M.gguf"
set SAMPLE_DIR=%CD%\sample_design
set OUTPUT_DIR=%CD%\sample_design\design_matrices
set API_URL=http://127.0.0.1:1234/v1/chat/completions
set MODEL_ID=gemma-2-9b-it

:: Create output directory
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

:: Check if LM Studio CLI is available
if not exist %LM_STUDIO_PATH% (
    echo Error: LM Studio CLI not found at %LM_STUDIO_PATH%
    echo Please edit the script to set the correct path
    exit /b 1
)

:: Start the LM Studio server
echo Starting LM Studio server...
%LM_STUDIO_PATH% server start

:: Load the model
echo Loading model...
%LM_STUDIO_PATH% load %MODEL_PATH%

:: Wait for server to start
echo Waiting for server to start...
timeout /t 5 /nobreak

:: Process each sample info file
for %%f in ("%SAMPLE_DIR%\*_sample_info.txt") do (
    echo Processing file: %%f
    
    :: Extract GSE ID from filename
    set filename=%%~nf
    set gse_id=!filename:_sample_info=!
    echo GSE ID: !gse_id!
    
    :: Create prompt file
    echo Analyze these samples and assign them to groups based on their experimental conditions. > "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    echo IMPORTANT: Use ONLY standardized group names as follows: >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    echo - If there is only one control and treatment such as 'WT' and 'KD', use the terms 'control' and 'treatment' respectively. >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    echo - For treatment samples that are not repliacates ^( i.e. more than one type of treatment^): treatment1, treatment2, etc. >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    echo For example, if there are 2 control samples and 3 treatment samples, they should be labeled as: >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    echo control, control, treatment1, treatment1, treatment2, treatment2, treatment3, treatment3 >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    echo If you are unsure, label them as unknown1, unknown2, etc. >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    echo These standardized names are REQUIRED for batch correction with ComBat-seq. >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    echo Please provide a tab-separated list with two columns: Sample_geo_accession and Group. >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    echo Sample information: >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    
    :: Extract both GSM IDs and sample titles - Fixed delimiter
    for /f "tokens=1,2 delims=	" %%a in ('type "%%f" ^| findstr /v "Sample_geo_accession"') do (
        echo GSM ID: %%a, Sample Title: %%b >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    )
    
    :: Include overall design information if available
    set design_file="%SAMPLE_DIR%\!gse_id!_overall_design.txt"
    if exist !design_file! (
        echo. >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
        echo Overall design information: >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
        type !design_file! >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    )
    
    echo. >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    
    echo The output must include the column headers 'Sample_geo_accession' and 'Group'. >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    echo Use ONLY the standardized group names ^(control1, treatment1, etc.^) as described above. >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    echo ONLY REPLY WITH THE TAB-SEPARATED LIST OF SAMPLE_GEO_ACCESSION AND GROUP IN ENGLISH ONLY, WITHOUT ANY OTHER TEXT OR LANGUAGES. >> "%OUTPUT_DIR%\!gse_id!_prompt.txt"
    
    :: Create a temporary file for the prompt content
    set "prompt_content="
    for /f "usebackq delims=" %%a in ("%OUTPUT_DIR%\!gse_id!_prompt.txt") do (
        set "prompt_content=!prompt_content!%%a\n"
    )
    
    :: Create JSON payload for the API request
    echo { > "%OUTPUT_DIR%\!gse_id!_payload.json"
    echo   "model": "%MODEL_ID%", >> "%OUTPUT_DIR%\!gse_id!_payload.json"
    echo   "messages": [ >> "%OUTPUT_DIR%\!gse_id!_payload.json"
    echo     { >> "%OUTPUT_DIR%\!gse_id!_payload.json"
    echo       "role": "user", >> "%OUTPUT_DIR%\!gse_id!_payload.json"
    echo       "content": "!prompt_content!" >> "%OUTPUT_DIR%\!gse_id!_payload.json"
    echo     } >> "%OUTPUT_DIR%\!gse_id!_payload.json"
    echo   ], >> "%OUTPUT_DIR%\!gse_id!_payload.json"
    echo   "temperature": 0.3, >> "%OUTPUT_DIR%\!gse_id!_payload.json"
    echo   "max_tokens": 3000 >> "%OUTPUT_DIR%\!gse_id!_payload.json"
    echo } >> "%OUTPUT_DIR%\!gse_id!_payload.json"
    
    :: Run curl to communicate with the LM Studio API
    echo Running API request for !gse_id!...
    curl -s -X POST "%API_URL%" ^
         -H "Content-Type: application/json" ^
         -d "@%OUTPUT_DIR%\!gse_id!_payload.json" ^
         -o "%OUTPUT_DIR%\!gse_id!_response.json"
    
    :: Extract the content from the JSON response using PowerShell
    echo Extracting content from response...
    powershell -Command "$response = Get-Content '%OUTPUT_DIR%\!gse_id!_response.json' -Raw | ConvertFrom-Json; $response.choices[0].message.content | Out-File -Encoding utf8 '%OUTPUT_DIR%\!gse_id!_group_assignments.txt'"
    
    :: Copy to design matrix file
    if exist "%OUTPUT_DIR%\!gse_id!_group_assignments.txt" (
        copy "%OUTPUT_DIR%\!gse_id!_group_assignments.txt" "%OUTPUT_DIR%\!gse_id!_design_matrix.txt" >nul
        echo Created design matrix for !gse_id!
    ) else (
        echo Failed to generate design matrix for !gse_id!
    )
)

:: Stop the LM Studio server
echo Stopping LM Studio server...
%LM_STUDIO_PATH% server stop

echo Processing complete!
endlocal 