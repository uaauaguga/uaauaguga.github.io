---
layout: post
title:  "Notes on SMTP"
date:   2021-07-31 21:30:12 +0800
usemathjax: true
categories: jekyll update
---

## Play with SMTP

### 用telnet和SMTP服务器交互在命令行发送邮件

```bash
# connect to 163's SMTP server
telnet smtp.163.com 25
## Trying 220.181.12.18...
## Connected to smtp.163.com.
## Escape character is '^]'.
## 220 163.com Anti-spam GT for Coremail System (163com[20141201])
HELO fs # a shaking is required prior to sending email 
## 250 OK
auth login # login SMTP server
## 334 dXNlcm5hbWU6 # this is a base64 encoding, validate this in python with base64.b64decode("dXNlcm5hbWU6"), which return b'username:'
xxxxxxxxxxxx # based 64 encoded user name
## 334 UGFzc3dvcmQ6 # b'Password:'
xxxxxxxxxxxx # based 64 encoded user password
## 235 Authentication successful
mail from:<jinyf1998@163.com> # input source email address
## 250 Mail OK
rcpt to:<jyf20@mails.tsinghua.edu.cn> # input target email adress
## 250 Mail OK
data # input email content
## 354 End data with <CR><LF>.<CR><LF>
subject: hello hawaii
from: <jinyf1998@163.com>
to: <jyf20@mails.tsinghua.edu.cn>

tell me, is some thing eluding you?
sunshine?
Is this not what you expected to see?
If you wanna find out what is behind this cold eyes ...
.
## 250 Mail OK queued as smtp14,EsCowADnavwKVwVhQPiWxQ--.36143S2 1627740013
# 这时jyf20@mails.tsinghua.edu.cn会收到来自jinyf1998@163.com的邮件
```

### 在python中利用SMTP发送邮件

- Send a email contains plain text

```python
# Import required libraries
import email
import smtplib
from email.utils import formataddr
from email.mime.text import MIMEText

# define content of the email
message = MIMEText('What is behind this cold eyes', 'plain', 'utf-8')
message['From'] = formataddr(("sender",sender))
message['To'] =  formataddr(("receiver",receiver))
message['Subject'] = 'Test Python SMTP'

# user information for login
sender = 'xxx'
receiver = 'xxx'
password = 'xxx'

smtpObj = smtplib.SMTP() 
# connect to SMTP server
smtpObj.connect("smtp.163.com", 25)    
smtpObj.login(sender,password)  

# send the email
smtpObj.sendmail(sender, receiver, message.as_string())
```

- Attatch a file to the email

```python
from email.utils import formataddr
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email import encoders

message = MIMEMultipart()
message['From'] = formataddr(("sender",sender))
message['To'] = formataddr(("receiver",receiver))
message['Subject'] = 'Test Python SMTP with attachment'

message.attach(MIMEText('hey you', 'plain', 'utf-8'))

with open('lyric.txt.gz', 'rb') as f:
    mime = MIMEBase('application', 'x-gzip', filename='lyric.txt.gz')
    mime.add_header('Content-Disposition', 'attachment', filename='lyric.txt.gz')
    mime.set_payload(f.read())
    encoders.encode_base64(mime)
    message.attach(mime)

# ... send the message as for plain text
```