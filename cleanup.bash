#/bin/bash

remove_cb()
{
  awk '{gsub(/_CB/,""); print}' $1 > _$1
  mv _$1 $1
}

remove_cd()
{
  awk '{gsub(/_CD/,""); print}' $1 > _$1
  mv _$1 $1
}

remove_params_global_bq()
{
  echo $1  
  awk '{gsub(/PARAMS_GLOBAL_BQ/,"PARAMS_GLOBAL"); print}' $1 > _$1
  mv _$1 $1
}


for i in `ls *_bq.f90`; do 
  remove_cb $i
  remove_cd $i
  remove_params_global_bq $i
  echo $i
done

rm -f *_cb.f90 *~
rm -f params_global_bq.f90

