#!/bin/bash

glslc particle.frag -o ../cmake-build-debug/shaders/frag_p.spv
glslc particle.vert -o ../cmake-build-debug/shaders/vert_p.spv

glslc comp.comp -o ../cmake-build-debug/shaders/comp.spv
glslc comp_it.comp -o ../cmake-build-debug/shaders/comp_it.spv
glslc vertex_shader.vert -o ../cmake-build-debug/shaders/vert.spv
glslc fragment_shader.frag -o ../cmake-build-debug/shaders/frag.spv