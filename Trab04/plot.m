set yrange [:] reverse
plot "VELOCIDADES.bin" binary format="%float" array=601x401 scan=yx with image
