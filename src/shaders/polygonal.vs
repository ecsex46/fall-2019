/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#version 330 core


#pragma once
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform mat4 normalMat;

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;

out vec3 vNormal;
out vec3 vPos;


void main()
{   
    vNormal = (normalMat * vec4(normal, 0.0)).xyz;
	vPos = (model * vec4(position, 1.f)).xyz;
    gl_Position = projection * view * model * vec4(position, 1.0);   
}

