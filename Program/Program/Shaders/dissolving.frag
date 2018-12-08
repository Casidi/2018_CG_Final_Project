#version 400 core

in vec3 fragNormal;
in vec2 TexCoord;
in vec3 fragPos;

out vec4 FragColor;

uniform sampler2D myTexture;
uniform sampler2D noiseTexture;

uniform float dissolvingThresh;

uniform vec3 lightPos;
uniform vec3 viewPos;

void main() {
	vec3 lightColor = vec3(1.0, 1.0, 1.0);

	vec3 norm = normalize(fragNormal);
	vec3 lightDir = normalize(lightPos - fragPos);

	float ambientStrength = 0.2;
	vec3 ambient = ambientStrength * lightColor;

	float diffuseStrenth = 0.8;
	float diff = max(dot(norm, lightDir), 0.0);
	vec3 diffuse = diffuseStrenth * diff * lightColor;

	float specularStrength = 0.5;
	vec3 viewDir = normalize(viewPos - fragPos);
	vec3 reflectDir = reflect(-lightDir, norm);
	float spec = pow(max(dot(viewDir, reflectDir), 0.0), 100);
	vec3 specular = specularStrength * spec * lightColor;

	float distance = length(lightPos - fragPos);
	float attenuation = 1.0 / (1.0 + distance * distance);

	vec3 result = ambient + diffuse + specular;
	vec3 resultTextured = ((ambient + diffuse)*vec3(texture(myTexture, TexCoord)) + specular) * attenuation;
	
	float noise = (texture(noiseTexture, TexCoord).x
					+ texture(noiseTexture, TexCoord).y
					+ texture(noiseTexture, TexCoord).z) / 3.0;
	if (noise < dissolvingThresh)
		discard;
	else if (noise < dissolvingThresh + 0.03)
		FragColor = vec4(1.0, 1.0, 0.0, 1.0);
	else
		FragColor = texture(myTexture, TexCoord);

	//FragColor = vec4(fragNormal, 1.0f);
	//FragColor = texture(myTexture, TexCoord);
	//FragColor = texture(noiseTexture, TexCoord);
	//FragColor = vec4(result, 1.0);
	//FragColor = vec4(resultTextured, 1.0);
}