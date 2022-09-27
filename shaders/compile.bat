::New and improved shader compilation batch script

for %%F in (*.vert, *.comp, *.frag) do (
    ::echo %%~nF.spv
    glslc %%F -o ../cmake-build-debug/shaders/%%~nF.spv
    glslc %%F -o %%~nF.spv
)