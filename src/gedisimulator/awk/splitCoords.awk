BEGIN{
  numb=0;
  n=-1;
}

($0&&($1!="#")){
  if((numb%maxPer)==0)n++;
  namen=sprintf("%s.%d.coords",root,n);
  print $0 >> namen

  numb++;
}

