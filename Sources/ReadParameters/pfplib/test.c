
void MAIN()
{
  int nx;
  char in[40];
  float array[1000];

  ifetch("nx",&nx);
  sfetch("in",in);

/*  rite(outfd,array,1000*sizeof(float));
*/
  printf("nx=%d\n",nx);
}
