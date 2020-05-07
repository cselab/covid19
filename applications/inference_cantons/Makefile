D = data/cantons/

default: run

run: samples.png tiles.png

samples.png: $D/state.pickle
	python3 -m korali.plotter --dir $D/_korali_samples --output $@ || true

tiles.png: $D/state.pickle plot.py
	./plot.py --tiles

$D/state.pickle: main.py ode.py data.py model.py
	./main.py

clean:
	rm -rf data
	rm -f samples.png prediction*.png data.{png,mp4} fit.{png,mp4} all.png

plot_map:
	./plot.py --map_data --map_fit --map_error --image --movie


.PHONY: default clean run plot_map
