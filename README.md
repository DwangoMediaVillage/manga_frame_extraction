manga-frame-extraction (MFE)
======================

## Build
```
cd MangaFrameExtraction
cmake ./
make
```

## Usage
###  Recursive Directory Scan Mode
```
./MFE -r -d /input/dir/ -o /output/dir/
```
###  Single Image Scan Mode
```
./MFE -s -f /input/file.jpg -o /output/dir/
```
###  Option
* set the maximum width of initialization before split execution.
```
-w 1000
```
* set the maximum height of initialization before split execution.
```
-h 1000
```
### Help
```
./MFE -h
```

## Require
*  OpenCV (https://github.com/Itseez/opencv)
*  Boost (http://www.boost.org 1.40 or later)

## License
MIT License, see [LICENSE](./LICENSE).

## Reference
1. 石井, 河村, 渡辺: “コミック画像のコマ分割処理における制御パラメータに関する検討”, 電子情報通信学会パターン認識・メディア理解研究会 PRMU2009-34 (IE2009-43, MI2009-34), pp.187-192, May 2009
