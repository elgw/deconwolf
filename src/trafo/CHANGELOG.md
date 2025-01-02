## 0.1.5

- Faster predictions: aborts directly when a single class has > 50% of
  the votes. Works in parallel regardless of row or column order of
  the input table.

## 0.1.4

- Added a minimal program to test that the library is installed and
  can be used in `test/minimal_example.c`.

- Hiding internal symbols in `libtrafo.so`, making it more tidy:

  ```
  $ nm libtrafo.so | grep ' T '
  0000000000004d00 T trafo_fit
  00000000000046d0 T trafo_free
  00000000000059f0 T trafo_importance
  00000000000056f0 T trafo_load
  00000000000049a0 T trafo_predict
  0000000000004750 T trafo_print
  0000000000005530 T trafo_save
  0000000000005480 T trafo_ut
  ```

## 0.1.3

- Added a few proxies for feature importance via `trafo_importance`.

- Fixed memory initialization to 0 after `realloc` (when needed).

- `trafo_cli` prints confusion matrices when the number of features is low.

## 0.1.1

- Added Entropy as a splitting criterion.
