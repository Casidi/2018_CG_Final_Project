#version 400 core

in vec3 fragNormal;
in vec2 TexCoord;
in vec3 fragPos;
in vec3 fragPosLocal;

out vec4 FragColor;

uniform sampler2D myTexture;
uniform vec3 lightPos;
uniform vec3 viewPos;
uniform float time;

void main() {
	vec3 norm = normalize(fragNormal);
	vec3 viewDir = normalize(viewPos - fragPos);

	FragColor = vec4(fragNormal, 1.0f);
	//FragColor = texture(myTexture, TexCoord);
	//FragColor = vec4(vec3(gl_FragCoord.z), 1.0);
}