---
layout: post
title:  "Notes on distributed file system"
date:   2022-08-30 19:52:08 +0800
usemathjax: true
categories: jekyll update
---

# Lustre


## Features

- See <https://lustre.ornl.gov/lustre101-courses/content/C1/L1/LustreIntro.pdf>

- Lustre server: servers for parallel storage
  - Object Storage Servers (OSS): Handles I/O requests for file data 
  - Object Storage Targets (OST): Block device used by OSS to store file data. Each OSS usually serves multiple OSTs
  - Metadata Storage Servers (MDS): Manages filenames and directories, file stripe locations, locking, ACLs, etc
  - Metadata Storage Targets (MDT):  Block device used by MDS to store metadata information 
  - MGS – Management server. Stores configuration information for one or more Lustre file systems.
  - MGT - Block device used by MGS for data storage 
- Lustre clients: Compute Nodes that access data from lustre server


## Resources

- https://wiki.lustre.org/Main_Page
- https://lustre.ornl.gov/lustre101-courses/