public class colorNote {
  
  float intensity;
  float mintensity;
  float frequency;
  color col;
  float[] spec;
  
  public colorNote()
  {
    intensity = 0;
    mintensity = 0;
    frequency = 0;
    col = color(0);
    spec = new float[10];
  }
  
  public colorNote(float i, float m, float f, color c, float[] s) 
  {
    intensity = i;
    mintensity = m;
    frequency = f;
    col = c;
    spec = s;
  }
}
