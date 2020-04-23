D = data/cantons/

default: run

run: samples.png prediction.png

samples.png: $D/state.pickle
	python3 -m korali.plotter --dir $D/_korali_samples --output $@ || true

prediction.png: $D/state.pickle
	./plot.py
	cp -v $D/figures/$@ $@

$D/state.pickle: main.py
	./main.py

clean:
	rm -rf data
	rm -f samples.png prediction.png

.PHONY: default clean run
