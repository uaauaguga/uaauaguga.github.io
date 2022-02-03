---
layout: post
title:  "Notes on ffmpeg"
date:   2022-02-03 10:06:32 +0800
usemathjax: true
categories: jekyll update
---

- ffmpeg documentation: <https://www.ffmpeg.org/documentation.html>

- Several tutorials
  - <https://ostechnix.com/20-ffmpeg-commands-beginners/>
  - <https://www.ruanyifeng.com/blog/2020/01/ffmpeg.html>
  - <https://www.zl-asica.com/2020/ffmpeg/>

### Crop image size in a video

```bash
ffmpeg -i input.mp4  -filter:v "crop=640:480:200:150" -y -v error output.mp4
# -filter:v mean a filter is applied to the video stream. crop=640:480:200:150 means crop a rect of size w:h:x:y
# if -y specified, overwrite existed file
# -v control the verbosity. 
# can be "quiet", "panic", "fatal", "error", "warning", "info", "verbose", "debug", "trace"
```

### Probe format and realted information of a video

```bash
ffprobe input.mp4 -v info -hide_banner
# -hide_banner: not print version and configuration of ffprobe
```

- Show duration:

```bash
ffprobe -v error -show_entries format=duration -of default=noprint_wrappers=1:nokey=1 input.mp4
# -show_entries entry_list  show a set of specified entries
# -of format alias for -print_format format, 
# set the output printing format default, compact, csv, flat, ini, json, xml
```

- Show number of stream

```bash
ffprobe -show_entries format=nb_streams -v 0 -of compact=p=0:nk=1 input.mp4
```
- Show bitrate:

```bash
# for video, 0 is index of the stream
ffprobe  -v error -select_streams v:0 -show_entries stream=bit_rate -of default=noprint_wrappers=1:nokey=1 input.mp4
# for audio, 0 is index of the stream
ffprobe  -v error -select_streams a:0 -show_entries stream=bit_rate -of default=noprint_wrappers=1:nokey=1 input.mp4
```

- A more complicated example
```bash
# https://superuser.com/questions/891665/ffprobe-show-entries-with-an-entry-name-that-uses-a-semicolon
ffprobe -v error -show_entries stream_tags=rotate:format=size,duration:stream=codec_name,bit_rate -of default=noprint_wrappers=1 input.mp4
# print as key value pair
```

### Slice a video

```bash
ffmpeg -ss 00:00:30.0 -i input.mp4 -c copy -t 00:00:10.0 output.mp4
# -ss: start time
# -t: duration (-to also works, which accepts stop time)
```



### Convert MP3 to MP4 by adding a static image

```bash
ffmpeg -loop 1 -i input.jpg -i input.mp3 -c:a copy -c:v libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuvj420p -y -shortest output.mp4
# -vf is an alias for -filter:v
# -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2": pad the image to make it dividable by 2
# -pix_fmt yuvj420p: some player can't parse the default pixcel format
```


### Concatenate multiple videos

```bash
for i in $(ls | grep '.mp3' );do echo "file $PWD/$i";done > filelist.txt
ffmpeg -f concat -safe 0 -i filelist.txt -c copy concatenated.mp3
```


### Extract images from video

```bash
mkdir -p  frames
ffmpeg -i input.mp4 -vf fps=1 frames/out-%03d.jpg
# fps=1: extract 1 image per second
```

### Extract audio from video

```bash
ffmpeg -i input.mp4 -vn output.mp3
```