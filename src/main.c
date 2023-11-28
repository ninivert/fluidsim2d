#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "solver.c"
#include "utils.h"

#include <glad/gl.h>
#include <GLFW/glfw3.h>

typedef struct {
  double x, y;
} Point;

struct {
  Point pos, pos_prev;
  bool pressed;
  int which_pressed;
} mouse;

static void error_callback(int error, const char* description) {
  fprintf(stderr, "Error: %s\n", description);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    glfwSetWindowShouldClose(window, GLFW_TRUE);
}

static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
  mouse.pressed = action == GLFW_PRESS;
  mouse.which_pressed = button;
}

static void mouse_position_callback(GLFWwindow* window, double x, double y) {
  mouse.pos_prev.x = mouse.pos.x;
  mouse.pos_prev.y = mouse.pos.y;
  mouse.pos.x = x;
  mouse.pos.y = y;
}

Point screen_to_world(Point screenpos, Point screenres, Point worldres) {
  return (Point) {
    linterp(screenpos.x, 0, screenres.x, 0, worldres.x),
    linterp(screenpos.y, 0, screenres.y, 0, worldres.y),
  };
}

int main(int argc, char** argv) {
  const uint nx = 200, ny = 100;
  const double source = 100.0, force = 5.0;

  // init simulation
  Sim* sim = sim_new(nx, ny, 0.1, 0.0, 0.0, 20);
  sim_print(sim);

  if (!glfwInit()) { fprintf(stderr, "failed to initialize glfw\n"); return -1; }
  glfwSetErrorCallback(error_callback);
  printf("initialized glfw\n");

  GLFWwindow* window = glfwCreateWindow(512, 512, "silly", NULL, NULL);
  if (!window) {
    fprintf(stderr, "failed to create glfw window\n");
    glfwTerminate();
    // TODO : de-alloc also the simulation
    return -1;
  }
  printf("initialized window\n");

  glfwMakeContextCurrent(window);
  if (!gladLoadGL(glfwGetProcAddress)) { fprintf(stderr, "failed to initialize opengl\n"); return -1; }
  printf("initialized opengl\n");
  
  glfwSetKeyCallback(window, key_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);
  glfwSetCursorPosCallback(window, mouse_position_callback);
  glfwSwapInterval(1);

  // initialize the texture we draw to
  GLuint txrid = 10; // GL texture ID
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, txrid);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 4);  // pixel data is word-aligned (every 4 bytes, i.e. 4*8=32bit to use floats)
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);  // no texture tiling
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);  // nearest filtering instead of linear
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_G, GL_RED);  // copy red component to green component
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_B, GL_RED);  // copy red component to blue component
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glDisable(GL_TEXTURE_2D);
  
  float* buffer = calloc(nx*ny, sizeof(float));
  float time = glfwGetTime(), timediff = 0.0;
  uint nframe = 0;

  while (!glfwWindowShouldClose(window)) {
    // update fps
    ++nframe;
    float newtime = glfwGetTime();
    timediff = newtime - time;
    time = newtime;

    // adjust for window resize		
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    glViewport(0, 0, width, height);

    // mouse interactions
    Point worldpos = screen_to_world(mouse.pos, (Point) {width, height}, (Point) {nx, ny});
    uint ix =      clampi((int) worldpos.x, 0, nx-1);
    uint iy = ny - clampi((int) worldpos.y, 0, ny-1) - 1;

    // reset dynamic source terms
    for (uint i = 0; i < nx*ny; ++i) { sim->vx_prev[i] = 0.0; }
    for (uint i = 0; i < nx*ny; ++i) { sim->vy_prev[i] = 0.0; }
    for (uint i = 0; i < nx*ny; ++i) { sim->rho_prev[i] = 0.0; }

    if (mouse.pressed) {
      if (mouse.which_pressed == GLFW_MOUSE_BUTTON_1) {
        // add a pointwise source at location of the mouse on the screen
        sim->rho_prev[IX(ix, iy)] = source;
      }
      if (mouse.which_pressed == GLFW_MOUSE_BUTTON_2) {
        // move the velocity field
        sim->vx_prev[IX(ix, iy)] = force * (mouse.pos.x - mouse.pos_prev.x);
        sim->vy_prev[IX(ix, iy)] = -force * (mouse.pos.y - mouse.pos_prev.y);
        // printf("vx=%f, vy=%f\n", A(sim->vx_prev, ix, iy), A(sim->vy_prev, ix, iy));
      }
    }

    // do computations
    sim_step(sim);

    // debug output
#ifdef DEBUG
    if (nframe % 60 == 0) {
      // printf("[nframe=%6d] screen (%.2f, %.2f), pixel (%d %d), rho=%f, vx=%f, vy=%f, source=%f\n", nframe, mouse.pos.x, mouse.pos.y, ix, iy, A(sim->rho, ix, iy), A(sim->vx, ix, iy), A(sim->vy, ix, iy), A(sim->source, ix, iy));
      printf("[nframe=%6d] frametime=%2.2fms (fps=%2.2f)\n", nframe, timediff*1e3f, 1.0/timediff);
    }
#endif
    
    // copy the simulation density to the texture buffer
    for (size_t i = 0; i < nx*ny; ++i) buffer[i] = (float) sim->rho[i];

    glClear(GL_COLOR_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // bind texture
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, txrid);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, nx, ny, 0, GL_RED, GL_FLOAT, buffer);
    // render single textured QUAD
    glColor3f(1.0,1.0,1.0);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0,0.0); glVertex2f(-1.0,-1.0);
    glTexCoord2f(1.0,0.0); glVertex2f(+1.0,-1.0);
    glTexCoord2f(1.0,1.0); glVertex2f(+1.0,+1.0);
    glTexCoord2f(0.0,1.0); glVertex2f(-1.0,+1.0);
    glEnd();
    // unbind texture
    glBindTexture(GL_TEXTURE_2D, 0);
    glDisable(GL_TEXTURE_2D);
    glFlush();

    glfwSwapBuffers(window);
    glfwPollEvents();
  }

  // free gl texture
  glDeleteTextures(1, &txrid);

  // free simulation resources
  sim_free(sim);
  free(buffer);

  // free glfw
  glfwDestroyWindow(window);
  glfwTerminate();
  
  return 0;
}

