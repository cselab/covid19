D = data/cantons/

default: run

run: samples.png prediction.png

samples.png: $D/state.pickle
	python3 -m korali.plotter --dir $D/_korali_samples --output $@ || true

prediction.png: $D/figures/prediction.png
	cp -v $< $@
	cp -v $D/figures/*.png .

$D/figures/prediction.png: main.py $D/state.pickle
	./plot.py

$D/state.pickle: main.py
	./main.py

clean:
	rm -rf data
	rm -f samples.png prediction.png

.PHONY: default clean run
