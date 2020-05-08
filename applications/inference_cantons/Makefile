D = data/cantons/

default: run

run: samples.png tiles.png

samples.png: $D/state.pickle
	python3 -m korali.plotter --dir $D/_korali_samples --output $@ || true

tiles.png: $D/state.pickle plot.py
	./plot.py --tiles --tiles_daily

$D/state.pickle: main.py ode.py data.py model.py
	./main.py

cleanplot:
	rm -vf samples.png prediction*.png {data,fit,error}*.{png,mp4}

cleandat:
	rm -rf data

clean: cleanplot cleandat

plot_map:
	./plot.py --map_data --map_fit --map_error --image --movie

plot_map_error:
	./plot.py --map_error --image

.PHONY: default clean run plot_map cleanplot cleandat plot_map_error
