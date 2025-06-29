module app;
import bindbc.glfw;
import bindbc.opengl;
import std.array : array;
import std.conv : to;
import std.csv;
import std.exception : enforce;
import std.file;
import std.math;
import std.random;
import std.stdio;

// # Grid resolution, should be square
int display = 1024;
// # Height of tip of pyramids (incorrect! the csv gives dist for steps of 0.1um)
// float top = 5;
// # Side length of area in which to plot, in microns
float length = 20;
// # Pyramids per square micron
float density = 0.6;

GLuint[2] programs;
GLuint[2] vaos;

void main() {
	enforce(loadGLFW() == GLFWSupport.glfw33, "Cant load GLFW_33");
	glfwSetErrorCallback(&glfw_error_callback);
	glfwInit();
	// char[] t = File("../tracing/heightCumulative.csv").byLine().front;
	// writeln(t);
	string csvText = readText("../tracing/heightCumulative.csv");
	float[] heightDist = array(csvReader!(float)(csvText).front);
	const uint pyramidCount = cast(uint) ceil(density * length * length);
	float[3][] peaks = genPeaks(heightDist, pyramidCount);

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	debug glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, true);
	GLFWwindow* window = glfwCreateWindow(display, display, "Pyramids", null, null);
	glfwMakeContextCurrent(window);
	if (window is null)
		writeln(glfwGetError(null));

	auto gl_version = loadOpenGL();
	enforce(gl_version == GLSupport.gl46, "Cant load OpenGL_46: " ~ gl_version.to!string);
	debug {
		glEnable(GL_DEBUG_OUTPUT);
		glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
		glDebugMessageCallback(&gl_error_callback, null);
		glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DEBUG_SEVERITY_NOTIFICATION, 0, null, false);
	}
	programs[0] = createShader(vertexPyramids, fragmentPyramids);
	programs[1] = createShader(vertexPlane, fragmentPlane);
	vaos[0] = createPyramid();
	vaos[1] = createPlane();

	setPeaks(peaks);
	setProjection();
	glfwSetScrollCallback(window, &scroll_callback);

	glClearColor(0, 0, 0, 1);
	glEnable(GL_DEPTH_TEST);
	while (!glfwWindowShouldClose(window)) {
		setCamera(window);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glUseProgram(programs[0]);
		glBindVertexArray(vaos[0]);
		glDrawElementsInstanced(GL_TRIANGLES, 12, GL_UNSIGNED_INT, cast(void*) 0, pyramidCount);

		glUseProgram(programs[1]);
		glBindVertexArray(vaos[1]);
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, cast(void*) 0);

		glfwSwapBuffers(window);
		glfwPollEvents();
		if (glfwGetKey(window, GLFW_KEY_ESCAPE))
			glfwSetWindowShouldClose(window, true);
	}
}

float[4][] perspectiveProjection(float aspectRatio, float horizontalFov = 2.1, // vertical fov 90°
	float nearplane = 0.01, float farplane = 200) {
	float a = 1.0 / tan(horizontalFov / 2.0);
	alias V = nearplane;
	alias A = farplane;
	alias s = aspectRatio;
	float z = -(A + V) / (A - V);
	float y = -(2.0 * A * V) / (A - V);
	return [[a, 0.0, 0.0, 0.0], [0.0, a * s, 0.0, 0.0], [0.0, 0.0, z, y], [0.0, 0.0, -1.0, 0.0]];
}

void setProjection() {
	foreach (program; programs) {
		glUseProgram(program);
		GLint loc = glGetUniformLocation(program, "projection");
		float[4][] proj = perspectiveProjection(1);
		glUniformMatrix4fv(loc, 1, true, proj[0].ptr);
	}
}

static float zoom = 1;
extern (C) void scroll_callback(GLFWwindow* window, double x, double y) nothrow {
	zoom += y * 0.1;
}

void setCamera(GLFWwindow* window) {
	double mx, my;
	glfwGetCursorPos(window, &mx, &my);
	mx /= display;
	my /= display;
	my = 1.5 - 3 * my;

	float radius = 2 * length / zoom;
	float z = cos(mx * 4 * PI);
	float x = sin(mx * 4 * PI);
	float y = sin(my * PI_4);
	float ny = cos(my * PI_4);

	float[4][] cam = [
		[z, y * -x, ny * x, radius * x], //
		[0, ny, y, radius * y], //
		[-x, y * -z, ny * z, radius * z], //
		[0, 0, 0, 1.0f], //
	];
	foreach (program; programs) {
		glUseProgram(program);
		GLint loc = glGetUniformLocation(program, "camera");
		glUniformMatrix4fv(loc, 1, true, cam[0].ptr);
	}
}

void setPeaks(float[3][] peaks) {
	GLuint peakbuffer;
	glBindVertexArray(vaos[0]);
	glCreateBuffers(1, &peakbuffer);
	glNamedBufferData(peakbuffer, peaks.length * 3 * float.sizeof, peaks[0].ptr, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, peakbuffer);
	glVertexAttribPointer(1u, 3, GL_FLOAT, false, 0, cast(void*) 0);
	glVertexAttribDivisor(1, 1);
	glEnableVertexAttribArray(1);
}

string vertexPlane = "#version 460
in vec3 pos;
out vec3 pos_world;
out vec4 pos_cam;
out vec3 normal;
out vec4 gl_Position;
uniform mat4 camera;
uniform mat4 projection;
void main(){
	normal = vec3(0,-1,0);
	pos_world = pos;
	pos_cam = inverse(camera) * vec4(pos_world, 1);
	gl_Position = projection * pos_cam;
}
	" ~ '\0';

string fragmentPlane = "#version 460
in vec3 normal;
in vec3 pos_world;
in vec4 pos_cam;
out vec4 color;
uniform mat4 camera;
uniform mat4 projection;
void main(){
	vec3 light = vec3(25, 25, 25);
	vec3 L = normalize(light - pos_cam.xyz/pos_cam.w);
	// vec3 light = normalize(vec3(1,1,1) - cam);
	float ambient = 0.3;
	float diffuse = max(0, dot(L, normal));
	color = (ambient + diffuse)*vec4(0.7,0,0.7,1);
}
	" ~ '\0';

string vertexPyramids = "#version 460
in vec3 pos;
in vec3 peak;
out vec3 pos_world;
out vec4 pos_cam;
out vec3 normal;
out vec4 gl_Position;
uniform mat4 camera;
uniform mat4 projection;
void main(){
	normal = normalize(pos);
	pos_world = (pos*peak.y) + vec3(peak.x, 0, peak.z);
	pos_cam = inverse(camera) * vec4(pos_world, 1);
	gl_Position = projection * pos_cam;
}
	" ~ '\0';
string fragmentPyramids = "#version 460
in vec3 normal;
in vec3 pos_world;
in vec4 pos_cam;
out vec4 color;
uniform mat4 camera;
uniform mat4 projection;
void main(){
	vec3 light = vec3(25, 25, 25);
	vec3 L = normalize(light - pos_cam.xyz/pos_cam.w);
	// vec3 light = normalize(vec3(1,1,1) - cam);
	float ambient = 0.3;
	float diffuse = max(0, dot(L, normal));
	color = (ambient + diffuse)*vec4(0.7,0.7,0,1);
}
	" ~ '\0';
GLuint createShader(string vertex, string fragment) {
	const(char*) vertexP = vertex.ptr;
	const(char*) fragmentP = fragment.ptr;

	writeln(vertex);
	writeln(fragment);

	GLuint vertShader = glCreateShader(GL_VERTEX_SHADER);
	GLuint fragShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(vertShader, 1, &vertexP, null);
	glShaderSource(fragShader, 1, &fragmentP, null);
	int succes;
	foreach (GLuint shader; [vertShader, fragShader]) {
		glCompileShader(shader);
		glGetShaderiv(shader, GL_COMPILE_STATUS, &succes);
		if (!succes) {
			int len;
			glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);
			char[] log = new char[len];
			glGetShaderInfoLog(shader, len, null, log.ptr);
			writeln(log);
			assert(0, log);
		}
	}
	GLuint shaderProgram = glCreateProgram();
	glAttachShader(shaderProgram, vertShader);
	glAttachShader(shaderProgram, fragShader);
	glLinkProgram(shaderProgram);
	glGetProgramiv(shaderProgram, GL_LINK_STATUS, &succes);
	if (!succes) {
		int len;
		glGetProgramiv(shaderProgram, GL_INFO_LOG_LENGTH, &len);
		char[] log = new char[len];
		glGetProgramInfoLog(shaderProgram, len, null, log.ptr);
		assert(0, log);
	}
	glDeleteShader(vertShader);
	glDeleteShader(fragShader);
	return shaderProgram;
}

GLuint createPyramid() {
	enum halfWidth = 1.0 / tan(54.7 / 180 * PI); // slope 54.7°?
	alias hw = halfWidth;
	float[] pos = [
		0, 1, 0, // peak
		-hw, 0, -hw, // bot left
		hw, 0, -hw, // bot right
		-hw, 0, hw, // top left
		hw, 0, hw, // top right
	];
	uint[] indices = [
		0, 1, 2, // bottom
		0, 3, 1, // left
		0, 4, 3, // top
		0, 2, 4, // right
	];
	return createVAO(pos, indices);
}

GLuint createPlane() {
	float w = length / 2.0f;
	float[] pos = [
		-w, 0, -w, // bot left
		w, 0, -w, //  bot right
		-w, 0, w, // top left
		w, 0, w, // top right
	];
	uint[] indices = [
		0, 1, 2, //
		2, 1, 3, //
	];
	return createVAO(pos, indices);
}

GLuint createVAO(float[] pos, uint[] indices) {
	GLuint vao, vboPos, ebo;
	glCreateVertexArrays(1, &vao);
	glCreateBuffers(1, &vboPos);
	glCreateBuffers(1, &ebo);

	glNamedBufferData(vboPos, pos.length * float.sizeof, pos.ptr, GL_STATIC_DRAW);
	glNamedBufferData(ebo, cast(long)(indices.length * uint.sizeof), indices.ptr, GL_STATIC_DRAW);

	glBindVertexArray(vao);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glBindBuffer(GL_ARRAY_BUFFER, vboPos);
	glVertexAttribPointer(0u, 3, GL_FLOAT, false, 0, cast(void*) 0);
	glEnableVertexAttribArray(0);

	glBindVertexArray(0);
	return vao;
}

float[3][] genPeaks(float[] dist, uint count) {
	float[3][] heights = new float[3][count];
	foreach (i; 0 .. count) {
		float x = uniform01() * length - length / 2;
		float z = uniform01() * length - length / 2;
		float yProb = uniform01();
		float y = sampleDist(dist, yProb);
		heights[i] = [x, y, z];
	}
	return heights;
}

float sampleDist(float[] dist, float zProb) {
	foreach_reverse (i, prob; dist) {
		if (dist[i] <= zProb) {
			return i / 10.0f; // um
		}
	}
	assert(0);
}

extern (C) void glfw_error_callback(int type, const char* description) nothrow {
	try {
		writefln("GLFW Exception %d: %s", type, description.to!string);
	} catch (Exception e) {
	}
}

extern (System) void gl_error_callback(GLenum source, GLenum type, GLuint errorID, GLenum severity,
	GLsizei length, const GLchar* message, const void* userParam) nothrow {
	try {
		const string[uint] sources = [
			GL_DEBUG_SOURCE_API: "OpenGL API", GL_DEBUG_SOURCE_WINDOW_SYSTEM: "Window System API",
			GL_DEBUG_SOURCE_SHADER_COMPILER: "Shader Compiler",
			GL_DEBUG_SOURCE_THIRD_PARTY: "Third Party",
			GL_DEBUG_SOURCE_APPLICATION: "Source Application", GL_DEBUG_SOURCE_OTHER: "Miscellaneous"
		];
		const string[uint] types = [
			GL_DEBUG_TYPE_ERROR: "Error ╮(. ❛ ᴗ ❛.)╭",
			GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: "Deprecated usage",
			GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR: "Undefined behaviour",
			GL_DEBUG_TYPE_PORTABILITY: "System portability",
			GL_DEBUG_TYPE_PERFORMANCE: "Performance Issues",
			GL_DEBUG_TYPE_MARKER: "\"Command stream annotation\"",
			GL_DEBUG_TYPE_PUSH_GROUP: "\"Group pushing\"",
			GL_DEBUG_TYPE_POP_GROUP: "\"Group popping\"", GL_DEBUG_TYPE_OTHER: "Miscellaneous"
		];
		const string[uint] severities = [
			GL_DEBUG_SEVERITY_HIGH: "High", GL_DEBUG_SEVERITY_MEDIUM: "Medium",
			GL_DEBUG_SEVERITY_LOW: "Low",
			GL_DEBUG_SEVERITY_NOTIFICATION: "Notification (Miscellaneous)"
		];

		writeln("Opengl Exception #", errorID);
		if (source !in sources || type !in types || severity !in severities)
			assert(0);
		writeln("\tSource: ", sources[source]);
		writeln("\tType: ", types[type]);
		writeln("\tSeverity: ", severities[severity]);
		writeln("\tMessage: ", message);
	} catch (Exception e) {
		try
			writeln(e);
		catch (Exception e) {
		}
	}
}
