# meradgen-cpp-final-dev

Working fork of `meradgen-cpp-final` for Approach A experiments
(sample at \(P=0\), freeze kinematics, `weight_at` for \(P=\pm1\)).

Do **not** treat this as the production library until deliberately
promoted. Production remains `meradgen-cpp-final/`.

See `NOTES.md` and `DESIGN_WEIGHT_AT.md`.

```bash
cmake -S meradgen-cpp-final-dev -B meradgen-cpp-final-dev/build \
  -DMERADGEN_BUILD_TESTS=ON
cmake --build meradgen-cpp-final-dev/build -j
./meradgen-cpp-final-dev/build/weight_at_smoke
```
