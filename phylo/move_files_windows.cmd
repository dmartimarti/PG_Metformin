FOR /f "tokens=* delims=" %a in (evo_strains.txt) DO move ".\assemblies\%a" ".\assemblies\evo_strains\"