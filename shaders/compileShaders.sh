#!/bin/bash
#New and improved shader compilation batch script

#for file in *.comp *.vert *.frag; do
#  glslc "$file" -o ../cmake-build-debug/shaders/"${file%.*}".spv
  #glslc "$file" -o "${file%.*}".spv
#done

for files in *.comp *.frag *.vert; do
  file=$(cat $files);
  extension="${files##*.}";

  while read line;
  do
    if [[ "$line" =~ "#include" ]]; then
      filename=${line/"#include "/""};
      incl="#include $filename";
      loaded=$(cat "$filename");
      file=${file/"$incl"/"$loaded"};
    fi
  done <<< "$file"

  echo "$file" >> tmp.$extension;

  glslc tmp.$extension -o ../cmake-build-debug/shaders/"${files%.*}".spv;
  rm tmp.$extension;
  #glslc "$file" -o "${file%.*}".spv
done