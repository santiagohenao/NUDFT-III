# NUDFT-III

Non Uniform Discrete Fourier Transform applied to OGLE-III Type 1 Cepheid photometry database.

## Usage

(linux)

`$ export JULIA_NUM_THREADS=$(nprocs)`

`$julia -p auto Main.jl n`

where n is the desired number of stars

## Results

Output file `results.dat`. 5 columns: star id, execution time, period, ephemeris and period quality factor (a.k.a. |F(T)|^2 on period T). Can be loaded on Julia with `read_float_table` from `FileOperations.jl`.

### Benckmark

With `n=1000` on Intel Core i3 legacy generation, Windows, quad core, Julia v1.1, CPU time was 805.928s and elapsed (calculation and output) time was 214.733s giving 375.3% of time improvement.

## Data References

* Udalski, Szymański, Soszyński and Poleski, 2008, _The Optical Gravitational Lensing Experiment. Final Reductions of the OGLE-III Data_, [Acta Astron., 58, 69](http://acta.astrouw.edu.pl/Vol58/n2/a_58_2_1.html).

* Soszyński et al., 2008a, _The Optical Gravitational Lensing Experiment. The OGLE-III Catalog of Variable Stars. I. Classical Cepheids in the Large Magellanic Cloud_, [Acta Astron., 58, 163](http://acta.astrouw.edu.pl/Vol58/n3/a_58_3_2.html)