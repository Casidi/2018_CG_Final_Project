#version 400 core

in vec3 fragNormal;
in vec2 TexCoord;
in vec3 fragPos;
in vec3 fragPosLocal;
in vec3 fragMatColor;

out vec4 FragColor;

uniform sampler2D myTexture;
uniform vec3 lightPos;
uniform vec3 viewPos;
uniform float time;

void main() {
	vec3 lightColor = vec3(1.0, 1.0, 1.0);

	vec3 norm = normalize(fragNormal);
	vec3 lightDir = normalize(lightPos - fragPos);

	float ambientStrength = 0.1;
	vec3 ambient = ambientStrength * lightColor;

	float diffuseStrenth = 0.5;
	float diff = max(dot(norm, lightDir), 0.0);
	vec3 diffuse = diffuseStrenth * diff * lightColor;

	float specularStrength = 0.5;
	vec3 viewDir = normalize(viewPos - fragPos);
	vec3 reflectDir = reflect(-lightDir, norm);
	float spec = pow(max(dot(viewDir, reflectDir), 0.0), 10);
	vec3 specular = specularStrength * spec * lightColor;

	float distance = length(lightPos - fragPos);
	float attenuation = 1.0 / pow(distance, 1);

	vec3 resultTextured = ambient * vec3(texture(myTexture, TexCoord))
							+ diffuse * vec3(texture(myTexture, TexCoord)) * attenuation 
							+ specular * attenuation;
	vec3 resultNoTexture = ambient + diffuse * fragMatColor * attenuation + specular * attenuation;

	//FragColor = vec4(fragNormal, 1.0f);
	//FragColor = vec4(fragMatColor, 1.0f);
	//FragColor = texture(myTexture, TexCoord);
	//FragColor = vec4(vec3(gl_FragCoord.z), 1.0);
	//FragColor = vec4(resultTextured, 1.0f);
	FragColor = vec4(resultNoTexture, 1.0f);
}