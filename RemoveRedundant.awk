BEGIN{i=0;}

{
if(NR>1 && $1!=b[1])
  {
  k=0;
  for(j=0;j<i;j++)
    {
    if(psig[j]+tsig[j] <= psig[k]+tsig[k] && flag[j]<=flag[k])
      k=j;
    if(flag[j]<flag[k])
      k=j;
    }
  print(line[k]);
  i=0;
  };

a=$0;
split(a,b);
line[i]=a;
mech[i]=b[12];
flag[i]=b[13];
psig[i]=b[14];
tsig[i]=b[15];
i++;
}

END{print(a)}
