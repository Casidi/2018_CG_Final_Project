#version 400 core

in vec3 fragNormal;
in vec2 TexCoord;
in vec3 fragPos;
in vec3 fragPosLocal;

out vec4 FragColor;

uniform sampler2D myTexture;
uniform sampler2D depthTexture;
uniform vec3 lightPos;
uniform vec3 viewPos;
uniform float time;

float linearize(float a) {
	float n = 0.001;                                // the near plane
	float f = 100.0;                               // the far plane
	return (2.0 * n) / (f + n - a * (f - n));  // convert to linear values       
}

float checkIntersect(float a) {
	a = linearize(a);
	float currentDepth = linearize(gl_FragCoord.z);
	if(a > currentDepth) {
		if (a - currentDepth < 0.01)
			return 1.0 - 100 * (a - currentDepth);
		else
			return 0.0;
	} else {
		return 0.0;
	}
}

void main() {
	vec3 norm = normalize(fragNormal);
	vec3 viewDir = normalize(viewPos - fragPos);
	
	float raw = texture2D(depthTexture, gl_FragCoord.xy / 512.0).r;
	if (gl_FragCoord.z > raw)
		discard;

	float Pole = pow((fragPosLocal.y + 1.0)/2, 2);
	float Rim = 1.0 - abs(dot(norm, viewDir));
	float Intersect = checkIntersect(raw);	

	vec4 _Color = vec4(16, 59, 159, 4) / 255.0;
    float glow = max(max(Intersect, Rim), Pole);
	vec4 glowColor = mix(_Color, vec4(1.0), glow);

	vec2 mainUV = vec2(TexCoord.x * 5, TexCoord.y * 9);
	vec4 mainTex = texture2D(myTexture, mainUV);
	float multiplier = (sin(fragPosLocal.y*6 + time*4) + 1.0) / 2.0;
	vec4 pulsedTexColor  = multiplier * mainTex * _Color;

	vec4 outColor = pulsedTexColor + glowColor;

	//FragColor = vec4(fragNormal, 1.0f);
	//FragColor = vec4(Rim);
	FragColor = outColor;
	//FragColor = vec4(vec3(linearize(texture(depthTexture, gl_FragCoord.xy / 512.0).r)), 1.0);
	//FragColor = texture(myTexture, TexCoord);
	//FragColor = vec4(vec3(gl_FragCoord.z), 1.0);

	//if(raw > gl_FragCoord.z && raw - gl_FragCoord.z < 0.00003)
	//	FragColor = vec4(1,0,0, 1);
}