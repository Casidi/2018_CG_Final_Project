#version 400 core

layout (location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aTexCoord;

out vec3 fragNormal;
out vec2 TexCoord;
out vec3 fragPos;

uniform mat4 modelMatrix;
uniform mat4 viewMatrix;
uniform mat4 projectionMatrix;

void main() {
	gl_Position = projectionMatrix * viewMatrix * modelMatrix * vec4(aPos, 1.0);
	fragNormal = mat3(transpose(inverse(modelMatrix))) * aNormal;
	TexCoord = aTexCoord;
	fragPos = vec3(modelMatrix * vec4(aPos, 1.0));
}