/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#version 330 core


#define PI 3.14159265359
#define MAX_NUM_EMITTER_TRIANGLES 40 // Max number of emitter triangles allowed (tuned for A4)
uniform float emitterVertices[MAX_NUM_EMITTER_TRIANGLES*3*3]; // Need to specify a size at compile time (max size is 512 floats)

uniform int nbTriangles; // Use nbTriangles to index the buffer correctly
uniform vec3 lightIrradiance;
uniform vec3 albedo;
uniform vec2 windowSize; // [width, height]

uniform sampler2D cvTerm; // Creates a 2D texture we can sample from (range x,y = [0,1])

in vec3 vNormal;
in vec3 vPos;

out vec3 color;

// Compute edge (v1--v2) contribution
float getEdgeContrib(vec3 v1, vec3 v2, vec3 pos) {
	// Adapt your getEdgeContrib code from the offline part
	float value = 0.f;
// TODO(A4): Implement this
	return value;
}


void main()
{	
	// 1) Extract vertices of triangles from `emitterVertices` buffer using `nbTriangles`
	// 2) Calculate G term
	// 3) Subtract modification term for G after extracting it from texture (use built-in `texture()` function)
	//	    e.g. `vec3 delta = texture(cvTerm, coords).xyz;`

	color = vec3(0);

    // TODO(A4): Implement this
}

