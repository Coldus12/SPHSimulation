::New and improved shader compilation batch script
::for %%F in (*.vert, *.comp, *.frag) do (
    ::echo %%~nF.spv
    ::glslc %%F -o ../cmake-build-debug/shaders/%%~nF.spv
    ::glslc %%F -o %%~nF.spv
::)

::Using multiline batch 'hack' from: https://stackoverflow.com/questions/7041069/passing-around-multi-line-strings
@echo off
setlocal EnableDelayedExpansion

::ls as in line separator
set ls=^


::Powershell script in batch (Puts text from #include-d file into the shader itself, and then compiles it)
set multiline=!ls!^
function parseInclude([string]$filename) {                                        !ls!^
                                                                                  !ls!^
  $file = Get-Content .\$filename;                                                !ls!^
  $raw_file = Get-Content .\$filename -Raw;                                       !ls!^
                                                                                  !ls!^
  foreach($line in $file) {                                                       !ls!^
    if ($line -match '#include') {                                                !ls!^
      $incl_filename = $line -replace \"#include \",\"\";                         !ls!^
      $incl_file = Get-Content .\$incl_filename -Raw;                             !ls!^
      $raw_file = $raw_file -replace \"#include $incl_filename\",\"$incl_file\";  !ls!^
    }                                                                             !ls!^
  }                                                                               !ls!^
                                                                                  !ls!^
  $extension = (Get-Item $filename ).Extension;                                   !ls!^
  $out_name = (Get-Item $filename ).Basename;                                     !ls!^
                                                                                  !ls!^
  echo \"$raw_file\" ^| Out-File -FilePath .\temp$extension -Encoding ascii;      !ls!^
  $t = Get-Content temp$extension -Raw;                                           !ls!^
  glslc .\temp$extension -o \"../cmake-build-debug/shaders/$out_name.spv\";       !ls!^
  #glslc .\temp$extension -o \"$out_name.spv\";                                   !ls!^
                                                                                  !ls!^
  Remove-Item .\temp$extension;                                                   !ls!^
  #Remove-Item \"$out_name.spv\";                                                 !ls!^
}

set ml2=!ls!^
function testFunc([string]$param1) {!ls!^
    echo \"ez powershellbol van: $param1\";!ls!^
}

:: Loops through all shader files while calling the powershell function on them
for %%F in (*.vert, *.comp, *.frag) do (
    call :func %%F
)
::@powershell !multiline!parseInclude(\"comp.comp\");
exit 0

:: Batch function
:func
@powershell !multiline!parseInclude(\"%~1\");

::Sources:
::  https://superuser.com/questions/1126265/windows-batch-powershell-command-inside-for-loop-not-working
::  https://stackoverflow.com/questions/127318/is-there-any-sed-like-utility-for-cmd-exe
::  https://stackoverflow.com/questions/7041069/passing-around-multi-line-strings
::  https://learn.microsoft.com/en-us/powershell/module/microsoft.powershell.core/about/about_regular_expressions?view=powershell-7.2
::  https://learn.microsoft.com/en-us/powershell/module/microsoft.powershell.core/about/about_functions?view=powershell-7.2
::  https://superuser.com/questions/1126265/windows-batch-powershell-command-inside-for-loop-not-working
::  https://stackoverflow.com/questions/12503871/removing-path-and-extension-from-filename-in-powershell