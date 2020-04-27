D = data/cantons/

default: run

run: samples.png prediction.png

samples.png: $D/state.pickle
	python3 -m korali.plotter --dir $D/_korali_samples --output $@ || true

prediction.png: $D/state.pickle plot.py
	./plot.py

$D/state.pickle: main.py ode.py data.py
	./main.py

clean:
	rm -rf data
	rm -f samples.png prediction*.png data.{png,mp4} fit.{png,mp4} all.png

.PHONY: default clean run
