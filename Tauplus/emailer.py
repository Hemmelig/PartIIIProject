#!/usr/bin/env python3

import smtplib

from email.mime.text import MIMEText

def sendEmail(address, message):

    msg = MIMEText(message)

    me = "conor.bacon@gmail.com"
    you = address
    
    msg['Subject'] = 'Code is done'
    msg['From'] = me
    msg['To'] = you

    server = smtplib.SMTP('localhost')
    server.sendmail(me, [you], msg.as_string())
    server.quit()
