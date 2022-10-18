import smtplib
from datetime import datetime
import os 

# def email_users(): 
gmail_user = 'seti.lofar@gmail.com'
gmail_password = 'ezturltallafrdmj'

time_stamp = datetime.now()
ip = 'Swedish Station'

sent_from = gmail_user
to = ['ojohnson@tcd.ie']
subject = 'TurboSETI has finished runnng'
body = ("TurboSETI search finished at %s on the %s" % (time_stamp, ip))

email_text = """\
From: %s
To: %s
Subject: %s

%s
""" % (sent_from, ", ".join(to), subject, body)

try:
    server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
    server.ehlo()
    server.login(gmail_user, gmail_password)
    server.sendmail(sent_from, to, email_text)
    server.close()

    print('Email sent!')
except:
    print('Something went wrong...')