mkdir -p data
mkdir -p figures

# Normal Mode I
#PYTHONPATH=../../../build:$PYTHONPATH python3 make_synthetic_with_deaths.py \
#    -D 3.7 -Z 3.7 -Zl 0.423 -Y 2.2 -alpha 0.3 -mu 0.4


# Normal Mode II
#PYTHONPATH=../../../build:$PYTHONPATH python3 make_synthetic_with_deaths.py \
#    -D 5.2 -Z 5.2 -Zl 0.423 -Y 2.2 -alpha 0.3 -mu 0.4 -r0sir 2.25 -r0seir 4.0 -r0seiir 7.0 -r0saphire 1.25 -r0seiru 2.25


# Upper Bounds
PYTHONPATH=../../../build:$PYTHONPATH python3 make_synthetic_with_deaths.py \
    -D 9.0 -Z 5.2 -Zl 0.423 -Y 3.0 -alpha 0.3 -mu 0.4 -r0sir 3.25 -r0seir 4.5 -r0seiir 8.5 -r0saphire 1.50 -r0seiru 3.25
