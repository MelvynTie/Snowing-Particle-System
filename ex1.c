//https://gamedevelopment.tutsplus.com/tutorials/create-a-cozy-snowy-night-scene-using-particle-effects--gamedev-3120
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#ifdef MACOSX
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif
#include "SOIL.h"

/******************************************************************************
* Constant
******************************************************************************/

#define RENDER 0
#define DEFAULT_LAST_TIME -1
#define ENABLER_COUNTER TRUE
#define MAX_NUM_PARTICLE 10000000
#define NANO_SECONDS 1000000000
#define DEG_TO_RAD 0.017453293
#define PI_ 3.142
#define RUN_SPEED  10
#define BUFFER_STRING 200

#define GL_GPU_MEM_INFO_TOTAL_AVAILABLE_MEM_NVX 0x9048
#define GL_GPU_MEM_INFO_CURRENT_AVAILABLE_MEM_NVX 0x9049

/******************************************************************************
* Utility Definition
******************************************************************************/

typedef enum {
  FALSE,
  TRUE
} Boolean;

/******************************************************************************
* Particle Properties Constant
******************************************************************************/

size_t total_particle = 1000;
size_t total_alive = 0;
GLdouble particle_life = 4.0;
GLdouble velocity_offset = 0.0;
GLdouble emission_offset = 0.1f;
GLdouble emission_rate;
GLdouble gravity_offset = 0.0;
GLdouble wind_offset = 0.0;
GLdouble snow_range = 100;
GLdouble snow_r = 0.9, snow_g = 0.9, snow_b = 0.9, snow_a = 1.0;

GLfloat forceX = 0.0;
GLfloat forceY = 0.0;
GLfloat gravityX = 0.0;
GLfloat gravityY = 9.81;
GLfloat windX = 1.0;
GLfloat windY = 0.0;

GLfloat particle_mass = 10;

int render_method = 0;
int emitter_enabled = 1;
int wall_enabled = 1;
int overlap_enabled = 1;
int emit_enabled = 1;
int axisEnabled = 0;

size_t current_overlap = 0;
size_t max_overlap = 3000;

GLfloat snowflake_size = 1;

/******************************************************************************
* Particle definition
******************************************************************************/

typedef struct {
  GLfloat p[3];
  GLfloat v[3];
  GLfloat a[3];
  GLfloat c[4];

  GLfloat rotation;
  GLfloat rotationSpeed;

  GLfloat mass;
  GLfloat life;
  GLfloat size;
  Boolean alive;

} Particle;

Particle p[MAX_NUM_PARTICLE];

typedef struct {
  GLfloat p[3];
  Boolean overlap;
} Terrain;

Terrain t[MAX_NUM_PARTICLE];

/******************************************************************************
* Window Constant
******************************************************************************/

GLint width = 800, height = 600;

/******************************************************************************
* Camera Constant
******************************************************************************/

GLdouble mlat = 0.0;  /* Mouse look offset angles */
GLdouble mlon = 0.0;

GLfloat  eyex = 0.0;  /* Eye point*/
GLfloat  eyey = 0.0;
GLfloat  eyez = 700.0;

GLfloat  centerx = 0.0; /* View point*/
GLfloat  centery = 0.0;
GLfloat  centerz = 0.0;

GLfloat  upx = 0.0; /* Camera point*/
GLfloat  upy = -1.0;
GLfloat  upz = 0.0;

GLfloat  minX = -200.0;
GLfloat  maxX = 200.0;
GLfloat  minY = -200.0;
GLfloat  maxY = 200.0;
GLfloat  minZ = -200.0;
GLfloat  maxZ = 200.0;

int drag_x_origin;
int drag_y_origin;
int dragging = 0;

GLuint axisList;

int AXIS_SIZE = 100;


/******************************************************************************
* TEXT definition
******************************************************************************/

GLfloat TEXT_X_POSITION = 0.99;
GLfloat TEXT_Y_POSITION = 0.95;

char buffer_string[BUFFER_STRING];

/******************************************************************************
* Timer Constant
******************************************************************************/

Boolean counter_enabled = ENABLER_COUNTER;
long timePassedNanos = 0;
long lastTime = DEFAULT_LAST_TIME;
long accumulatedTimeNano = 0;

/******************************************************************************
* FPS Constant
******************************************************************************/

float fps = 60;
int frames = 0;
long firstFrameTime = 0;
int fpsRefreshTimeNanos = 500 * 1000 * 1000; //each 500ms

/******************************************************************************
* Prototypes
******************************************************************************/

void emitParticle();

void SnowingUpdater();

void renderer();

double getUniformRandom();

void displayMenu();

void counter();
long getCurrentNanoTime();
float getTimeAccumulatedSeconds();
float getTimePassedSeconds();

void showData();
void printString(void*,double,double,char*);

/******************************************************************************
* Counter
******************************************************************************/

void counter()
{
  if(!counter_enabled) return;

  if(lastTime == DEFAULT_LAST_TIME)
  {
    lastTime = getCurrentNanoTime();
    timePassedNanos = 0;
    fps = 0;
    frames = 0;
    firstFrameTime = lastTime;
  } else {
    long currentTime = getCurrentNanoTime();   //Get the current time
    timePassedNanos = currentTime - lastTime; //Time passed
    lastTime = currentTime;                   //Update last time
    accumulatedTimeNano += timePassedNanos;  //Accumulate time
    frames++;                                 //FPS
    long dt = currentTime - firstFrameTime;   //Calculate fps
    if(dt >= fpsRefreshTimeNanos) {
        fps = (float)(1000*frames)/(float)(dt / 1000000);
        frames = 0;
        firstFrameTime = currentTime;
    }
  }
}

/******************************************************************************
* Get Time Passed in Seconds
******************************************************************************/

float getTimePassedSeconds() { return (float) timePassedNanos / NANO_SECONDS;}

/******************************************************************************
* Get Accumulated Second Time
******************************************************************************/

float getTimeAccumulatedSeconds() { return (float) accumulatedTimeNano / NANO_SECONDS;}

/******************************************************************************
* Get Current Nano Time
******************************************************************************/

long getCurrentNanoTime()
{
  struct timespec ts;
  timespec_get(&ts, TIME_UTC);
  return (long)ts.tv_sec * 1000000000L + ts.tv_nsec;
}

/******************************************************************************
* Get Random number between range [0,1]
******************************************************************************/

double getUniformRandom() { return (rand()/(double)RAND_MAX);}

/******************************************************************************
* Init Terrain
******************************************************************************/

void initTerrain()
{
  for(size_t id=0; id<MAX_NUM_PARTICLE; id++)
  {
    t[id].p[0] = t[id].p[1] = t[id].p[2] = 0;
    t[id].overlap = FALSE;
  }
}

/******************************************************************************
* Emit particle
******************************************************************************/

// Keep emitting the particles
void emitParticle()
{
  // Find the delta time
  float dt = getTimePassedSeconds();
  // Control the emission_rate
  emission_rate = total_particle * emission_offset;
  // Max particles
  size_t max_new_particle = (size_t)(dt * emission_rate);
  size_t start_id = total_alive;
  size_t end_id = (start_id+max_new_particle < total_particle-1) ? (start_id+max_new_particle):(total_particle-1);

  // Generate the initial position of the snow
  size_t id;
  for(id=start_id; id<end_id; id++)
  {
    p[id].p[0] = getUniformRandom() * (snow_range - (-snow_range)) - snow_range;
    p[id].p[1] = minY;
    p[id].p[2] = getUniformRandom() * (snow_range - (-snow_range)) - snow_range;

    // Move slowly in one direction
    p[id].v[0] = getUniformRandom() * PI_ * 2.0;
    p[id].v[1] = getUniformRandom() * PI_ * 2.0;
    p[id].v[2] = 0;

    p[id].c[0] = snow_r;
    p[id].c[1] = snow_g;
    p[id].c[2] = snow_b;
    p[id].c[3] = snow_a;

    // Get random angle [0,360]
    p[id].rotation = 0;

    p[id].mass = particle_mass;
    // life [3.0,4.0]
    p[id].life = 3.0 + getUniformRandom();
    p[id].alive = FALSE;
    p[id].size = getUniformRandom() * (3.0 - (0.0)) + 0.0;
  }

  // Wake particle up
  for(id=start_id; id<end_id; id++)
  {
    if(total_alive<total_particle)
    {
      p[total_alive].alive = TRUE;

      Particle temp = p[id];
      p[id] = p[total_alive];
      p[total_alive] = temp;

      total_alive++;
    }
  }
}

/******************************************************************************
* Check Collision
******************************************************************************/

int checkInRange(GLfloat x, GLfloat y, GLfloat z) {
  int xInRange;
  int yInRange;
  int zInRange;

  if((minX < x)&&(maxX > x)) {xInRange = 1;}
    else {xInRange = 0;}

  if((minY < y)&&(maxY > y)) {yInRange = 1;}
    else {yInRange = 0;}

  if((minZ < z)&&(maxZ > z)) {zInRange = 1;}
    else {zInRange = 0;}

  return xInRange && yInRange && zInRange;
}

/******************************************************************************
* Snowing Updater
******************************************************************************/

void SnowingUpdater()
{
  // Get time passed
  GLfloat dt = getTimePassedSeconds();

  float counter = (double)accumulatedTimeNano / NANO_SECONDS;

  size_t end_id = total_alive;

  for(size_t id=1; id<end_id; id++)
  {
    // Sum all external forces
    // Update Acceleration
    forceX = (gravityX + windX) * p[id].mass;
    forceY = (gravityY + gravity_offset + windY) * p[id].mass;

    p[id].a[0] += forceX;
    p[id].a[1] += forceY;

    p[id].a[0] /= p[id].mass;
    p[id].a[1] /= p[id].mass;

    windX = pow(sin(counter) * 0.5 + 0.5, 4) * 10 - wind_offset;
    windY = sin(counter * 100);

    // Update new position
    p[id].p[0] += p[id].v[0] * dt;
    p[id].p[1] += p[id].v[1] * dt;
    p[id].p[2] += 0;

    // Update new velocity
    p[id].v[0] += (p[id].a[0] + velocity_offset) * dt;
    p[id].v[1] += (p[id].a[1] + velocity_offset) * dt;

    p[id].rotationSpeed += (p[id].a[0] / 10);

    p[id].rotation += p[id].rotationSpeed * dt;

    if(p[id].rotation >= 360.0) p[id].rotation = 0.0;

    p[id].life -= dt;

    int inRange = checkInRange(p[id].p[0], p[id].p[1], p[id].p[2]);
    if(!inRange)
    {
      // Overlapping Snow
      current_overlap = (current_overlap <= max_overlap) ? current_overlap : 0;

      if(t[current_overlap].overlap)
      {

        t[current_overlap].p[0] = p[id].p[0];
        t[current_overlap].p[1] = p[id].p[1];
        t[current_overlap].p[2] = p[id].p[2];

        current_overlap++;

      } else
      {
        t[current_overlap].p[0] = p[id].p[0];
        t[current_overlap].p[1] = p[id].p[1];
        t[current_overlap].p[2] = p[id].p[2];

        t[current_overlap].overlap = TRUE;
        current_overlap++;
      }


      // Kill particle
      if(total_alive>0)
      {
        p[id].p[0] = p[id].p[1] = p[id].p[2] = 0;
        p[id].v[0] = p[id].v[1] = p[id].v[2] = 0;
        p[id].a[0] = p[id].a[1] = p[id].a[2] = 0;
        p[id].c[0] = p[id].c[1] = p[id].c[2] = p[id].c[3] = 0.0;

        p[id].mass = particle_mass;
        p[id].life = 0;
        p[id].alive = FALSE;
        p[id].size = 0;

        Particle temp = p[id];
        p[id] = p[total_alive];
        p[total_alive] = temp;

        total_alive--;
        end_id = (total_alive < total_particle) ? total_alive : total_particle;
      }
    }
  }
}

/******************************************************************************
* Renderer
******************************************************************************/

void renderer()
{
  glClear(GL_COLOR_BUFFER_BIT);

  if(emitter_enabled)
  {
    // Draw Emitter
    glColor4f(0.8f, 0.8f, 0.8f, 0.8f);
    for(int i=-snow_range; i <= snow_range; i+= 2)
    {
      glBegin(GL_QUADS);
        glVertex3d(-snow_range, minY, snow_range);
        glVertex3d(snow_range, minY, snow_range);
        glVertex3d(snow_range, minY, -snow_range);
        glVertex3d(-snow_range, minY, -snow_range);
      glEnd();
    }
  }

  // Draw Point
  if(render_method == 0)
  {
    glEnable(GL_POINT_SMOOTH);
    glPointSize(2);
    for(size_t id=1; id<total_alive; id++)
    {
      if(p[id].alive)
      {
        glBegin(GL_POINTS);
        glColor4f(p[id].c[0],p[id].c[1],p[id].c[2],p[id].c[3]);
        glVertex3f(p[id].p[0], p[id].p[1], p[id].p[2]);
        glEnd();
      }
    }
  }
    // Draw Sprite
  else if(render_method == 1)
  {
    glPointSize(2);
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_POINT_SPRITE);

    // glColor4f(snow_r,snow_g,snow_b,snow_a);

    int bitmatArray = SOIL_load_OGL_texture("SnowParticle.png",SOIL_LOAD_RGBA,SOIL_CREATE_NEW_ID,SOIL_FLAG_MIPMAPS);
    glBindTexture(GL_TEXTURE_2D, bitmatArray);
     for(size_t id=1; id<total_alive; id++)
     {
       if(p[id].alive)
       {
         glBegin(GL_POINTS);
         glVertex3f(p[id].p[0], p[id].p[1], p[id].p[2]);
         glEnd();
       }
     }

     glDisable(GL_TEXTURE_2D);
     glDisable(GL_POINT_SPRITE);
  }
    // Draw texture
  else if(render_method == 2) {
    for(size_t id=1; id<total_alive; id++)
    {
      if(p[id].alive)
      {
        // Draw Triangle
        double cx = p[id].p[0];
        double cy = p[id].p[1] + p[id].size;

        double bx = cx - p[id].size / tan(60 * DEG_TO_RAD);
        double by = cy - p[id].size;

        double ax = cx + p[id].size / tan(60 * DEG_TO_RAD);
        double ay = cy - p[id].size;

        double tx = p[id].p[0];
        double ty = p[id].p[1] + p[id].size/2;

        double dx = tx;
        double dy = ty - p[id].size;

        double ex = tx + p[id].size / tan(60 * DEG_TO_RAD);
        double ey = ty;

        double fx = tx - p[id].size / tan(60 * DEG_TO_RAD);
        double fy = ty;

        glColor4f(p[id].c[0],p[id].c[1],p[id].c[2],p[id].c[3]);
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();

        glTranslatef(p[id].p[0], p[id].p[1], p[id].p[2]);
        glRotated(p[id].rotation,0.0,1.0,0.0);
        glTranslatef(-p[id].p[0], -p[id].p[1], -p[id].p[2]);

        glBegin(GL_TRIANGLES);
          glVertex3f(ax, ay, p[id].p[2]); // vertex 1
          glVertex3f(bx, by, p[id].p[2]); // vertex 2
          glVertex3f(cx, cy, p[id].p[2]); // vertex 3
        glEnd();

        glBegin(GL_TRIANGLES);
          glVertex3f(dx, dy, p[id].p[2]); // vertex 4
          glVertex3f(ex, ey, p[id].p[2]); // vertex 5
          glVertex3f(fx, fy, p[id].p[2]); // vertex 6
        glEnd();

        glPopMatrix();
      }
    }
  }

if(overlap_enabled)
{
  // Draw Terrain
  glEnable(GL_POINT_SMOOTH);
  glPointSize(5);
  for(size_t id=1; id<max_overlap; id++)
  {
    if(t[id].overlap)
    {
      glBegin(GL_POINTS);
      glVertex3f(t[id].p[0], t[id].p[1], t[id].p[2]);
      glEnd();
    }
  }
}


if(wall_enabled)
{
  // Draw Floor
  glColor4f(0.1f, 0.1f, 0.1f, 0.1f);
  for(int i=minX; i <= maxX; i+= 40)
  {
    glBegin(GL_LINES);
      glVertex3d(minX, maxY, i);
      glVertex3d(maxX, maxY, i);
    glEnd();

    glBegin(GL_LINES);
      glVertex3d(i, maxY, minZ);
      glVertex3d(i, maxY, maxZ);
    glEnd();
  }

  // Draw Left
  glColor4f(0.1f, 0.1f, 0.1f, 0.1f);
  for(int i=minZ; i <= maxZ; i+= 40)
  {
    glBegin(GL_LINES);
      glVertex3d(maxX, minY, i);
      glVertex3d(maxX, maxY, i);
    glEnd();

    glBegin(GL_LINES);
      glVertex3d(maxX, i, minZ);
      glVertex3d(maxX, i, maxZ);
    glEnd();
  }

  // Draw Right
  glColor4f(0.1f, 0.1f, 0.1f, 0.1f);
  for(int i=minZ; i <= maxZ; i+= 40)
  {
    glBegin(GL_LINES);
      glVertex3d(minX, minY, i);
      glVertex3d(minX, maxY, i);
    glEnd();

    glBegin(GL_LINES);
      glVertex3d(minX, i, minZ);
      glVertex3d(minX, i, maxZ);
    glEnd();
  }

  // Draw Ceiling
  glColor4f(0.1f, 0.1f, 0.1f, 0.1f);
  for(int i=minZ; i <= maxZ; i+= 40)
  {
    glBegin(GL_LINES);
      glVertex3d(minX, minY, i);
      glVertex3d(maxX, minY, i);
    glEnd();

    glBegin(GL_LINES);
      glVertex3d(i, minY, minZ);
      glVertex3d(i, minY, maxZ);
    glEnd();
  }

  // Draw Back
  glColor4f(0.1f, 0.1f, 0.1f, 0.1f);
  for(int i=minX; i <= maxX; i+= 40)
  {
    glBegin(GL_LINES);
      glVertex3d(i, minY, minZ);
      glVertex3d(i, maxY, minZ);
    glEnd();

    glBegin(GL_LINES);
      glVertex3d(minX, i, minZ);
      glVertex3d(maxX, i, minZ);
    glEnd();
  }
}

}

/******************************************************************************
* Main Function
******************************************************************************/

void display()
{
  // Load the same Matrix
   glLoadIdentity();

  gluLookAt(eyex, eyey, eyez,
            centerx, centery, centerz,
            upx, upy, upz);

  glRotated(mlon, 1.0, 0.0, 0.0);
  glRotated(mlat, 0.0, 1.0, 0.0);

  // Background Colour
  glClearColor(0.0, 0.0, 0.0, 1.0);
  // Clear the screen
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if(counter_enabled) counter();

  if(emit_enabled) emitParticle();

  SnowingUpdater();
  renderer();
  showData();

  if(axisEnabled) glCallList(axisList);

  glutPostRedisplay();
  glutSwapBuffers();
}

/******************************************************************************
* Make Axes
******************************************************************************/
void makeAxes() {
  axisList = glGenLists(1);
  glNewList(axisList, GL_COMPILE);
      glLineWidth(2.0);
      glBegin(GL_LINES);
      glColor3f(1.0, 0.0, 0.0);       // X axis - red
      glVertex3f(0.0, 0.0, 0.0);
      glVertex3f(AXIS_SIZE, 0.0, 0.0);
      glColor3f(0.0, 1.0, 0.0);       // Y axis - green
      glVertex3f(0.0, 0.0, 0.0);
      glVertex3f(0.0, AXIS_SIZE, 0.0);
      glColor3f(0.0, 0.0, 1.0);       // Z axis - blue
      glVertex3f(0.0, 0.0, 0.0);
      glVertex3f(0.0, 0.0, AXIS_SIZE);
    glEnd();
  glEndList();
}


/******************************************************************************
* Keyboard
******************************************************************************/

void keyboard(unsigned char key, int x, int y)
{
  switch (key) {
    case 27: // Escape Key
        exit(0);
      break;
    /* Emission Rate */
    case 'e':
      emission_offset -= 0.01;
      break;
    case 'E':
      emission_offset += 0.01;
      break;
    /* Emission Rate */

    /* Initial Velocity */
    case 'v':
      velocity_offset -= 0.1;
      break;
    case 'V':
      velocity_offset += 0.1;
      break;
    /* Initial Velocity */

    /* Life Time */
    case 't':
      particle_life -= 0.1;
      break;
    case 'T':
      particle_life += 0.1;
      break;
    /* Life Time */

    /* Particle Number */
    case 'P':
      if(total_particle != 0) total_particle -= 1000;
      break;
    case 'p':
      total_particle += 1000;
      break;
    /* Particle Number */

    /* Gravity */
    case 'g':
      gravity_offset -= 1;
      break;
    case 'G':
      gravity_offset += 1;
      break;
    /* Gravity */

    /* Wind */
    case 'w':
      wind_offset -= 1;
      break;
    case 'W':
      wind_offset += 1;
      break;
    /* Wind */
   }
  glutPostRedisplay();
}

/******************************************************************************
* Cursor Keyboard
******************************************************************************/

void cursor_keys(int key, int x, int y)
{
  switch (key) {

    case GLUT_KEY_UP:
      eyex = sin(mlon*DEG_TO_RAD) *  RUN_SPEED + eyex;
      eyez = cos(mlon*DEG_TO_RAD) *  RUN_SPEED + eyez;
      break;
    case GLUT_KEY_DOWN:
      eyex = eyex - sin(mlon*DEG_TO_RAD) *  RUN_SPEED;
      eyez = eyez - cos(mlon*DEG_TO_RAD) *  RUN_SPEED;
      break;
    case GLUT_KEY_LEFT:
      eyex = sin((mlon+90)*DEG_TO_RAD) *  RUN_SPEED + eyex;
      eyez = cos((mlon+90)*DEG_TO_RAD) *  RUN_SPEED + eyez;
      break;
    case GLUT_KEY_RIGHT:
      eyex = eyex - sin((mlon+90)*DEG_TO_RAD) * RUN_SPEED;
      eyez = eyez - cos((mlon+90)*DEG_TO_RAD) * RUN_SPEED;
      break;
  }
}

/******************************************************************************
* Mouse Click
******************************************************************************/

void mouse_click(int button, int state, int x, int y)
{
    if(button == GLUT_LEFT_BUTTON) {
        if(state == GLUT_DOWN) {
            dragging = 1;
            drag_x_origin = x;
            drag_y_origin = y;
        }
        else
            dragging = 0;
    }
}

/******************************************************************************
* Mouse Motion
******************************************************************************/

void mouse_motion(int x, int y)
{
    if(dragging) {
        mlon += (y - drag_y_origin)*0.3;
        mlat += (x - drag_x_origin)*0.3;
        drag_x_origin = x;
        drag_y_origin = y;
    }
}

/******************************************************************************
* Reshape
******************************************************************************/

void reshape(int w, int h)
{
  glClearColor(0.9, 0.9, 0.9, 1.0);
  glViewport(0, 0, (GLsizei)w, (GLsizei)h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45, (GLfloat)w / (GLfloat)h, -1, 50);
  glMatrixMode(GL_MODELVIEW);

  width = w;
  height = h;
}

/******************************************************************************
*  Show Data
******************************************************************************/

void showData()
{
  sprintf(buffer_string,"Particle number: %lu Alive: %lu  Emission Rate: %.2f FPS: %.2f ",total_particle, total_alive, emission_rate, fps);
  printString(GLUT_BITMAP_HELVETICA_12, TEXT_X_POSITION, TEXT_Y_POSITION,buffer_string);

  sprintf(buffer_string,"Gravity: -%.2f ms^-1 Wind: %.2f ms^-1",gravityY+gravity_offset, windX);
  printString(GLUT_BITMAP_HELVETICA_12, TEXT_X_POSITION, TEXT_Y_POSITION-0.05,buffer_string);

  float local_t = getTimeAccumulatedSeconds();
  sprintf(buffer_string,"Life: %.2f Time: %.2fs",particle_life, local_t);
  printString(GLUT_BITMAP_HELVETICA_12, TEXT_X_POSITION, TEXT_Y_POSITION-0.1,buffer_string);

  char* s = (char*) glGetString(GL_VENDOR);
  printString(GLUT_BITMAP_HELVETICA_12, TEXT_X_POSITION, TEXT_Y_POSITION-0.15,s);

  s = (char*) glGetString(GL_RENDERER);
  printString(GLUT_BITMAP_HELVETICA_12, TEXT_X_POSITION, TEXT_Y_POSITION-0.2,s);

}

/******************************************************************************
*  Print String
******************************************************************************/

// Performance bottleneck
void printString(void* font, double x, double y, char* string)
{
  glColor4f(1.0,1.0,1.0,1.0);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0.0, 0.0, 0.0, 0.0, -1.5, 1.5);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glRasterPos3f(-x,y,0.0);
  int len;
  len = (int)strlen(string);

  for(int i=0; i<len; i++) {
    glutBitmapCharacter(font, string[i]); }

  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);

}

/******************************************************************************
* Create Menu
******************************************************************************/

void menu (int menuentry)
{
  switch (menuentry) {
    case 1:
        accumulatedTimeNano = 0;
      break;
    case 2:
      eyex = 0.0;  /* Eye point*/
      eyey = 0.0;
      eyez = 700.0;

      centerx = 0.0; /* View point*/
      centery = 0.0;
      centerz = 0.0;

      upx = 0.0; /* Camera point*/
      upy = -1.0;
      upz = 0.0;

      break;
    case 3:
      eyex = 0.0;  /* Eye point*/
      eyey = 0.0;
      eyez = 1000.0;

      centerx = 0.0; /* View point*/
      centery = 0.0;
      centerz = 0.0;

      upx = 0.0; /* Camera point*/
      upy = -1.0;
      upz = 0.0;
      break;
    case 4:
      eyex = 0.0;  /* Eye point*/
      eyey = 0.0;
      eyez = 200.0;

      centerx = 0.0; /* View point*/
      centery = 100.0;
      centerz = 0.0;

      upx = 0.0; /* Camera point*/
      upy = -1.0;
      upz = 0.0;
    case 5:
      render_method = 0;
      break;
    case 6:
      render_method = 2;
      break;
    case 7:
      render_method = 1;
      break;
    case 8:
      axisEnabled = (axisEnabled) ? 0:1 ;
      break;
    case 9:
      emit_enabled = (emit_enabled) ? 0:1 ;
      break;
    case 10:
      emitter_enabled = (emitter_enabled) ? 0:1 ;
    break;
    case 11:
      wall_enabled = (wall_enabled) ? 0:1 ;
    break;
    case 12:
      overlap_enabled = (overlap_enabled) ? 0:1 ;
    break;
    case 13: exit(0);
  }
}

/******************************************************************************
* Main Function
******************************************************************************/

int main(int argc, char *argv[])
{
  srand(time(NULL));

  glutInit(&argc, argv);
  glutInitDisplayMode (GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
  glutInitWindowSize(width, height);
  glutInitWindowPosition(100, 100);
  glutCreateWindow("COMP37111 Particles");

  glutCreateMenu (menu);
  glutAddMenuEntry ("Rest Time", 1);
  glutAddMenuEntry ("Default View", 2);
  glutAddMenuEntry ("Distant View", 3);
  glutAddMenuEntry ("Close View", 4);
  glutAddMenuEntry ("Render(Points)", 5);
  glutAddMenuEntry ("Render(Triangles)", 6);
  glutAddMenuEntry ("Render(Image)", 7);
  glutAddMenuEntry ("On/Off Axis", 8);
  glutAddMenuEntry ("On/Off Emission", 9);
  glutAddMenuEntry ("On/Off Emitter", 10);
  glutAddMenuEntry ("On/Off Wall", 11);
  glutAddMenuEntry ("On/Off Overlap", 12);
  glutAddMenuEntry ("Exit", 13);
  glutAttachMenu (GLUT_RIGHT_BUTTON);

  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(cursor_keys);
  glutMotionFunc(mouse_motion);
  glutMouseFunc(mouse_click);
  glutReshapeFunc(reshape);

  makeAxes();
  initTerrain();

  glutMainLoop();
  exit(1);
}
