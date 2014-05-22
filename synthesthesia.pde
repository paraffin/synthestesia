import ddf.minim.*;
import ddf.minim.signals.*;
import ddf.minim.analysis.*;
import ddf.minim.effects.*;

Minim minim;
AudioPlayer jingle;
AudioInput lineIn;
FFT fft;
WaveformRenderer waveform;

PGraphics pgb, pgm, pgh;
PImage img;

int bands = 12;
int colors = 256;
int buffersize;
float sampleRate;
int delay = 2;

float[] bass, mids, high;
colorNote[] bass_col;
colorNote[] mids_col, high_col;

int bass_octaves, mids_octaves, high_octaves;
// # of octaves below/above 440Hz to start bass
int low_oct;

float rotation;


void setup()
{
  size(800, 600, P2D);
  background(0);
  fill(#ffff00);
  ellipseMode(CENTER);
  noFill();
  colorMode(HSB, colors, colors, colors, colors);
  
  minim = new Minim(this);
  
  //lineIn = minim.getLineIn();
  //fft = new FFT(lineIn.bufferSize(),lineIn.sampleRate());

  //jingle = minim.loadFile(selectInput("Please choose a music file"),4096);

  jingle = minim.loadFile("pater-noster.mp3", 4096);  

  jingle.play();

  fft = new FFT(jingle.bufferSize(), jingle.sampleRate());
  fft.noAverages();
  fft.window(FFT.HAMMING);

  waveform = new WaveformRenderer();
  jingle.addListener(waveform);

  bass_octaves = 2;
  mids_octaves = 3;
  high_octaves = 2;
  low_oct = -4;

  bass_col = new colorNote[bass_octaves*bands];
  mids_col = new colorNote[mids_octaves*bands];
  //high_col = new colorNote[high_octaves*bands];
  
  for(int i=0; i<mids_col.length; i++) { mids_col[i] = new colorNote(); }
  //for(int i=0; i<high_col.length; i++) { high_col[i] = new colorNote(); }

  bass = new float[bass_octaves*bands];
  mids = new float[mids_octaves*bands];
  high = new float[high_octaves*bands];

  float bass_max, mids_max, high_max;
}

synchronized void draw()
{
  
  //background(0);
  fft.forward(jingle.mix);
  //fft.forward(lineIn.mix);

  mids_col = processSpectrum(-1, 1, bands);
  bass_col = processBass(-3, -1, bands);
  mids_col = reprocess(mids_col);
  drawBassCircle(bass_col);
  drawMids();
  //waveform.drawCircle(color(colors,colors,colors,0));
  //drawColorSpectrum(bass_col, height);
  blend(g, 5, 5, width-9, height-9, 0, 0, width, height, BLEND);
  filter(ERODE);
  
  
 
  
  
  //blend(g, 0, 4, width+1, height-8, 0, 0, width, height, BLEND);
}

// Take FFT array fft, process into full colorNote array
// - Take only octave range; 440Hz = octave0 and n_bands
// - Go one note at a time, querying spectrum with frequency range
// - Each note split into n_bands
// - One color for each note
// - Normalized intensity for each note: (max - min)/max
// - Frequency is weighted average over a note
// - Hue based on this frequency
// - Saturation: either 
//   - 100% or 
//   - 1 - Normalized integral over note
// - Brightness: either
//   - 100% or
//   - 
// - Alpha:
//   - 100% or
//   - 
colorNote[] processSpectrum(int low_oct, int high_oct, int n_bands)
{
  colorNote[] rv = new colorNote[(high_oct - low_oct)*12];
  for(int i=0; i<rv.length; i++) { rv[i] = new colorNote(); }
  
  for(int oct = low_oct; oct < high_oct; oct++)
  {
    float power = 0;
    float freq = 0;
    float intensity = 0;
    
    for(int i=0; i<12; i++)
    {
      int index = (oct-low_oct)*12 + i;
      rv[index].spec = new float[n_bands];
      
      float sum1 = 0;
      float sum2 = 0;
      float maxi = 0;
      float mini = 0;
      float inten;
      
      for(int j=0; j<n_bands; j++)
      {
        power = ( ((float) oct*12.00) + ((float) i) - 0.5 + ( (float) j / ((float) n_bands) ) ) / ((float) 12.00);
        freq = (float) Math.pow(2, power)*440;
        
        inten = fft.getFreq(freq);
        
        sum1 += inten*freq;
        sum2 += inten;
        
        maxi = max(maxi, inten);
        mini = min(mini, inten);
        
        rv[index].spec[j] = inten;
      }
      
      float avg_f = sum1/sum2;
      float normInt = (maxi-mini)/maxi;
      
      float avg_p = (float) (12*Math.log(avg_f/440) / Math.log(2)) - 12*oct;
      
      int hew = Math.round(map(avg_p, -0.5, 12.5, 0, colors));
      int sat = Math.round(map(1 - sum2/(maxi*n_bands), 0, 1, 0.5, colors));
      //int sat = Math.round(map(normInt, 0, 1, 0.5, colors));
      int bri = sat;
      int alp = colors;
      
      rv[index].col = color(hew,sat,bri,alp);
      rv[index].intensity = maxi;
      rv[index].mintensity = sum2/bands;
      rv[index].frequency = avg_f;
    }
  }
  return rv;
}

colorNote[] processBass(int low_oct, int high_oct, int n_bands)
{
  colorNote[] rv = new colorNote[(high_oct - low_oct)*12];
  for(int i=0; i<rv.length; i++) { rv[i] = new colorNote(); }
  
  for(int oct = low_oct; oct < high_oct; oct++)
  {
    float power = 0;
    float freq = 0;
    float intensity = 0;
    
    for(int i=0; i<12; i++)
    {
      int index = (oct-low_oct)*12 + i;
      rv[index].spec = new float[n_bands];
      
      float sum1 = 0;
      float sum2 = 0;
      float maxi = 0;
      float mini = 0;
      float inten;
      
      for(int j=0; j<n_bands; j++)
      {
        power = ( ((float) oct*12.00) + ((float) i) - 0.5 + ( (float) j / ((float) n_bands) ) ) / ((float) 12.00);
        freq = (float) Math.pow(2, power)*440;
        
        inten = fft.getFreq(freq);
        
        sum1 += inten*freq;
        sum2 += inten;
        
        maxi = max(maxi, inten);
        mini = min(mini, inten);
        
        rv[index].spec[j] = inten;
      }
      
      float avg_f = sum1/sum2;
      float normInt = (maxi-mini)/maxi;
      
      float avg_p = (float) (12*Math.log(avg_f/440) / Math.log(2)) - 12*oct;
      
      int hew = Math.round(map(avg_p, -0.5, 12.5, 0, colors));
      int sat = colors;//Math.round(map(1 - sum2/(maxi*n_bands), 0, 1, 0.8, colors));
      //int sat = Math.round(map(normInt, 0, 1, 0.5, colors));
      int bri = sat;
      int alp = colors;
      
      rv[index].col = color(hew,sat,bri,alp);
      rv[index].intensity = maxi;
      rv[index].mintensity = sum2/bands;
      rv[index].frequency = avg_f;
    }
  }
  return rv;
}

// Try some post-processing stuff that needs a processed spec
// Notably, saturation and brightness modified here
colorNote[] reprocess(colorNote[] spec)
{
  colorNote[] tmp = spec;
  color tmpc;
  int normInt;
  
  for (int i=0; i<spec.length; i++) { tmp[i] = spec[i]; }
  for(int i=1; i<spec.length-1; i++)
  {    
    normInt = Math.round(map( (2*spec[i].intensity - spec[i-1].mintensity - spec[i+1].mintensity) / spec[i].intensity,
                   0, 1, 0, colors));
    tmpc = spec[i].col;
    tmp[i].col = color(hue(tmpc), normInt, normInt, alpha(tmpc));
  }
  return tmp;
}

void drawBassCircle(colorNote[] spec)
{
  float max_inten = 0;
  int max_index = 0;
  float[] weights = new float[12];
  
  for(int i=0; i<spec.length; i++)
  {
    weights[i%12] += spec[i].intensity;
  }
  float sum = 0;
  for(int i=0; i<12; i++)
  {
    if(weights[i] > max_inten)
    {
      max_index = i;
      max_inten = spec[i].intensity;
    }
    sum += spec[i].intensity;
  }
  int hew = Math.round(hue(spec[max_index].col));
  int sat = Math.round( map( (max_inten - sum/spec.length)/max_inten, 0, 1, .25*colors, colors) );
  int bri = sat;
  int alp = colors;
  color col = color(hew,sat,bri,alp);
  
  ellipseMode(CENTER);
  fill(col);
  stroke(col);
  ellipse(width/2, height/2, 2*max_inten, 2*max_inten);
}

// Draw arcs to make round spectrum, rotate each one by (the same)
// angle given by bass circle
void drawMids()
{
  rotation += map(max(bass), 0, 300, PI/512, PI/128);
  float dTheta = TWO_PI/(bands*mids_col.length);
  ellipseMode(CENTER);
  
  for(int i=0; i<mids_col.length; i++)
  {
    fill(mids_col[i].col);
    stroke(mids_col[i].col);
    
    for(int j=0; j<bands; j++)
    {
      float rad = mids_col[i].spec[j]*15;
      arc(width/2, height/2, rad, rad, ((i*bands + j)*dTheta + rotation)%TWO_PI, ((i*bands+j+1)*dTheta + rotation)%TWO_PI);
    }
  }
  //blend(pgm, 1, 1, width-1, height-1, 0, 0, width, height, REPLACE);
} 
void drawHigh()
{
  //rotation += map(max(bass), 0, 300, PI/512, PI/32);
  float dTheta = TWO_PI/high_col.length;
  
  for(int i=0; i<high_col.length; i++)
  {
    float rad = high_col[i].intensity*15;
    fill(high_col[i].col);
    stroke(high_col[i].col);
    ellipseMode(CENTER);
    arc(width/2, height/2, rad, rad, (i*dTheta + rotation)%TWO_PI, ((i+1)*dTheta + rotation)%TWO_PI) ;
  }
  //blend(pgm, 1, 1, width-1, height-1, 0, 0, width, height, REPLACE);
} 

void drawBassCircleNoBuf()
{
  // find mean frequency in bass, set as color, draw circle to bass buffer
  float sum1 = 0;
  float sum2 = 0;
  for(int i=0; i<bass_col.length; i++)
  {
   sum1 += bass_col[i].intensity*i; 
   sum2 += bass_col[i].intensity;
  }
  int mean = Math.round(sum1/sum2);
  color col = bass_col[mean].col;
  
  ellipseMode(CENTER);
  fill(col);
  stroke(col);
  ellipse(width/2, height/2, 2*bass_col[mean].intensity, 2*bass_col[mean].intensity);
}
  
  



  // Draw given array of colors as vertical lines across the screen
  // with given height y
  synchronized void drawColorSpectrum(colorNote[] spec, int y)
  {
    for (int i=0; i<width/bands; i++)
    {
      int index = Math.round( (float)i*bands*(spec.length-1) / ((float) width) );
      stroke(spec[index].col);
      for(int j=0; j<bands; j++)
      {
        line(i*bands+j, y, i*bands+j, y-5*spec[index].spec[j]);
      }
    }
  }

  // Draw given array's values as vertical lines across the screen
  // with given height y
  synchronized void drawSpectrum(float[] spec, int y)
  {
    for (int i=0; i<width; i++)
    {
      int index = Math.round( (float)i*(spec.length-1) / ((float) width) );
      stroke(y, 256, 256);
      line(i, y, i, y-spec[index]);
    }
  }

colorNote[] fftToColors(int lo, int ox, colorNote[] prev_col)
{
  colorNote[] rv = new colorNote[ox*bands];
  
  float power;
  float freq_hi;
  float freq_lo;
  float intensity;
  
  for(int i=lo; i<lo+ox; i++)
  {
    for(int j=0; j<bands; j++)
    {
      int index = (i-lo)*bands+j;
      power = (((float) i*bands) + ((float) j) - 0.5) / ((float) bands);
      freq_lo = (float) Math.pow(2, power)*440;
      power = (((float) i*bands) + ((float) j) + 0.5) / ((float) bands);
      freq_hi = (float) Math.pow(2, power)*440;
      
      float sum1 = 0;
      float sum2 = 0;
      float max_in = 0;
      float tmpf;
      int count = 0;
      for(float k=freq_lo; k<freq_hi; k+=1.00)
      {
        tmpf = fft.getFreq(k);
        sum1 += tmpf*k;
        sum2 += tmpf;
        max_in = max(max_in, tmpf);
        count++;
      }
      float mean_f = sum1/sum2;
      
      intensity = (log(mean_f/440+1)+0.5)*max_in;
      
      // Get the float/mean equivalent of j so we can convert to hue
      float mean_j = (float) (bands*Math.log(mean_f/440) / Math.log(2)) - i*bands;
      
      float hew = map(mean_j, 0, bands, 0, colors);
      float sat = colors;//map(max_in * count / sum2, 0, 10, colors/1.5, colors);
      /*
      float bri = map(intensity, 0, 100, colors/4, 2*colors)
                  - 0.2*brightness(prev_col[(i-lo)*bands+j]);
      float alp = map(2*bri, colors/2, colors, 0, colors)+0.1*alpha(prev_col[(i-lo)*bands+j]);      
      */
      float bri = map(intensity, 0, 200, colors/4, colors);
      float alp = map(intensity, 0, 200, 0, colors);
      
      float rvi = 0.75*intensity + 0.25*prev_col[index].intensity;
      rv[(i-lo)*bands+j] = new colorNote(rvi, 0, mean_f, color(hew,bri,sat,alp), new float[1]);
    }
  }
  return rv;
}
