#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import run_subprocess
import selenium.common.exceptions
from rauth import OAuth1Session
import multiprocessing
import subprocess
import getpass
import sys
import os
import re
# import chromedriver_autoinstaller
import chromedriver_binary
from cryptography.fernet import Fernet
# selenium automation
from selenium import webdriver
# from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options

"""
Script to test access to authenticated resources via REST interface.
Written by Keith Jolley
Copyright (c) 2017, University of Oxford
E-mail: keith.jolley@zoo.ox.ac.uk

This file is part of Bacterial Isolate Genome Sequence Database (BIGSdb).

BIGSdb is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BIGSdb is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The test databases can be reached at https://pubmlst.org/test/.
To use these, sign up for a PubMLST account (https://pubmlst.org/site_accounts.shtml)
and link this account with the pubmlst_test_seqdef and pubmlst_test_isolates 
databases (https://pubmlst.org/site_accounts.shtml#registering_with_databases)
"""


'modified by adamkoziol'


class REST(object):

    def main(self):
        """
        Run the appropriate methods in the correct order
        """
        self.secret_finder()
        self.parse_access_token()
        self.get_session_token()
        self.parse_session_token()
        self.get_route()
        self.download_profile()
        self.find_loci()
        self.download_loci()

    def secret_finder(self):
        """
        Parses the supplied secret.txt file for the consumer key and secrets
        """
        secretlist = list()
        if os.path.isfile(self.secret_file):
            # Open the file, and put the contents into a list
            with open(self.secret_file, 'r') as secret:
                for line in secret:
                    secretlist.append(line.rstrip())
            # Extract the key and secret from the list
            self.consumer_key = secretlist[0]
            self.consumer_secret = secretlist[1]
        else:
            print('"Cannot find the secret.txt file required for authorization. '
                  'Please ensure that this file exists, and that the supplied consumer key is on the '
                  'first line, and the consumer secret is on the second line. '
                  'Contact keith.jolley@zoo.ox.ac.uk for an account, and the necessary keys')
            quit()

    def parse_access_token(self):
        """
        Extract the secret and token values from the access_token file
        """
        access_file = os.path.join(self.file_path, 'access_token')
        # Ensure that the access_token file exists
        if os.path.isfile(access_file):
            # Initialise a list to store the secret and token
            access_list = list()
            with open(access_file, 'r') as access_token:
                for line in access_token:
                    value, data = line.split('=')
                    access_list.append(data.rstrip())
            # Set the variables appropriately
            self.access_secret = access_list[0]
            self.access_token = access_list[1]
        else:
            print('Missing access_token')
            self.get_request_token()
            self.get_access_token()

    def create_encryption_key_file(self):
        # key generation
        key = Fernet.generate_key()

        # string the key in a file
        with open(self.credentials_key, 'wb') as file_key:
            file_key.write(key)

    def read_encryption_key(self):
        # Open the key file
        with open(self.credentials_key, 'rb') as file_key:
            key = file_key.read()
        # Use the generated key
        fernet = Fernet(key)
        return fernet

    def encrypt_credentials(self):

        login_user = input("Please enter your PubMLST username:\n")
        login_pass = getpass.getpass(prompt='Please enter your PubMLST password:\n').encode('utf-8').decode()
        # Write the credentials to file
        with open(self.credentials_file, 'w') as credentials:
            credentials.write(f'{login_user}\n{login_pass}')
        # Create an encryption key to encrypt the file
        self.create_encryption_key_file()
        fernet = self.read_encryption_key()
        # Open the original file to encrypt
        with open(self.credentials_file, 'rb') as file:
            original = file.read()
        # Encrypt the file
        encrypted = fernet.encrypt(original)
        # Open the file in write mode and write the encrypted data
        with open(self.credentials_file, 'wb') as encrypted_file:
            encrypted_file.write(encrypted)
        return login_user, login_pass

    def decrypt_credentials(self):

        # Read in and decrypt the encryption key
        fernet = self.read_encryption_key()
        # Open the encrypted file
        with open(self.credentials_file, 'rb') as enc_file:
            encrypted = enc_file.read()
        # Decrypt the file
        decrypted = fernet.decrypt(encrypted).decode()
        # Set the username and password from the decrypted, decoded string
        login_user, login_pass = decrypted.split('\n')
        return login_user, login_pass

    def get_access_token(self):
        """
        Request an access token
        """
        print("Authorizing access")
        # Set URL to use for the verification
        authorize_url = self.test_web_url + '&page=authorizeClient&oauth_token=' + self.request_token
        if not os.path.isfile(self.credentials_file):
            login_user, login_pass = self.encrypt_credentials()
        else:
            login_user, login_pass = self.decrypt_credentials()
            if not login_user or login_pass:
                login_user, login_pass = self.encrypt_credentials()
        try:
            verifier = self.selenium_authorise(authorize_url=authorize_url,
                                               login_user=login_user,
                                               login_pass=login_pass)
        except selenium.common.exceptions.SessionNotCreatedException:
            try:
                # Attempt to install the correct version of chromedriver
                run_subprocess('python -m pip install chromedriver-autoinstaller --trusted-host pypi.org '
                               '--trusted-host files.pythonhosted.org --trusted-host pypi.python.org ')
                import chromedriver_autoinstaller
                chromedriver_autoinstaller.install()
                verifier = self.selenium_authorise(authorize_url=authorize_url,
                                                   login_user=login_user,
                                                   login_pass=login_pass)
            except (ModuleNotFoundError, selenium.common.exceptions.SessionNotCreatedException):
                print('Visit this URL in your browser: ' + authorize_url)
                # Use the user input to set the verifier code
                verifier = input('Enter oauth_verifier from browser: ')
        # Create a new session
        session_request = OAuth1Session(consumer_key=self.consumer_key,
                                        consumer_secret=self.consumer_secret,
                                        access_token=self.request_token,
                                        access_token_secret=self.request_secret)
        # Perform a GET request with the appropriate keys and tokens
        r = session_request.get(self.access_token_url,
                                verify=False,
                                params={
                                    'oauth_verifier': verifier
                                })
        # If the status code is '200' (OK), proceed
        if r.status_code == 200:
            # Save the JSON-decoded token secret and token
            self.access_token = r.json()['oauth_token']
            self.access_secret = r.json()['oauth_token_secret']
            # Write the token and secret to file
            self.write_token('access_token', self.access_token, self.access_secret)

    def selenium_authorise(self, authorize_url, login_user, login_pass):

        # Setup chrome options
        chrome_options = Options()
        chrome_options.add_argument("--headless")  # Ensure GUI is off
        chrome_options.add_argument("--no-sandbox")
        chrome_options.add_argument("ignore-certificate-errors")
        # Start up browser
        browser = webdriver.Chrome(options=chrome_options)

        # Get page
        browser.get(authorize_url)
        assert "Log in" in browser.title

        # Find and fill inputs
        user = browser.find_element_by_name("user")
        user.clear()
        user.send_keys(login_user)
        password = browser.find_element_by_name("password_field")
        password.clear()
        password.send_keys(login_pass)
        # If the cookie compliance element is present, it will intercept the login click
        try:
            cookie_dismiss = browser.find_element_by_class_name('cc-compliance')
            cookie_dismiss.click()
        except (selenium.common.exceptions.ElementNotInteractableException,
                selenium.common.exceptions.NoSuchElementException,
                selenium.common.exceptions.InvalidElementStateException):
            pass
        login_button = browser.find_element_by_name("submit")
        login_button.click()
        authorize_button = browser.find_element_by_name("submit")
        authorize_button.click()

        assert "Authorize third-party client" in browser.title
        code = browser.find_element_by_xpath("/html/body/div[2]/div[2]/div/div[2]/p[3]/b").get_attribute("innerHTML")
        verifier = code.split()[2]
        browser.quit()
        print("Selenium Authorization Successful")
        return verifier

    def get_request_token(self):
        """
        Obtain a request token
        """
        print('Obtaining request token')
        try:
            os.remove(os.path.join(self.file_path, 'request_token'))
        except FileNotFoundError:
            pass
        # Create a new session
        session = OAuth1Session(consumer_key=self.consumer_key,
                                consumer_secret=self.consumer_secret)
        # Use the test URL in the GET request
        r = session.request(method='GET',
                            url=self.request_token_url,
                            verify=False,
                            params={'oauth_callback': 'oob'})
        # If the status code is '200' (OK), proceed
        if r.status_code == 200:
            # Save the JSON-decoded token secret and token
            self.request_token = r.json()['oauth_token']
            self.request_secret = r.json()['oauth_token_secret']
            # Write the token and secret to file
            self.write_token('request_token', self.request_token, self.request_secret)

    def get_session_token(self, attempt=0):
        """
        Use the accession token to request a new session token
        """
        # self.logging.info('Getting session token')
        # Rather than testing any previous session tokens to see if they are still valid, simply delete old tokens in
        # preparation of the creation of new ones
        try:
            os.remove(os.path.join(self.file_path, 'session_token'))
        except FileNotFoundError:
            pass
        # Create a new session
        session_request = OAuth1Session(self.consumer_key,
                                        self.consumer_secret,
                                        access_token=self.access_token,
                                        access_token_secret=self.access_secret)
        # Perform a GET request with the appropriate keys and tokens
        r = session_request.get(self.session_token_url,
                                verify=False)
        # If the status code is '200' (OK), proceed
        if r.status_code == 200:
            # Save the JSON-decoded token secret and token
            self.session_token = r.json()['oauth_token']
            self.session_secret = r.json()['oauth_token_secret']
            # Write the token and secret to file
            self.write_token('session_token', self.session_token, self.session_secret)
        # Any other status than 200 is considered a failure
        else:
            print('Failed:')
            print(r.json()['message'])
            if 'Invalid access token.  Generate new access token' in str(r.json()['message']) and not attempt:
                self.overwrite_access_token()

    def overwrite_access_token(self):
        self.get_request_token()
        self.get_access_token()
        self.get_session_token(attempt=1)

    def write_token(self, token_type, token, secret):
        """
        Write a token to file. Format is secret='secret'\n,token='token'
        :param token_type: The type of token. Options are 'request', 'session', and 'access'
        :param token: The string of the token extracted from the GET request
        :param secret:
        """
        # Open the file, and write the token and secret strings appropriately
        with open(os.path.join(self.file_path, token_type), 'w') as token_file:
            token_file.write('secret=' + secret + '\n')
            token_file.write('token=' + token + '\n')

    def parse_session_token(self):
        """
        Extract the session secret and token strings from the session token file
        """
        session_file = os.path.join(self.file_path, 'session_token')
        # Only try to extract the strings if the file exists
        if os.path.isfile(session_file):
            # Create a list to store the data from the file
            session_list = list()
            with open(session_file, 'r') as session_token:
                for line in session_token:
                    # Split the description e.g. secret= from the line
                    value, data = line.split('=')
                    # Add each string to the list
                    session_list.append(data.rstrip())
            # Extract the appropriate variable from the list
            self.session_secret = session_list[0]
            self.session_token = session_list[1]

    def get_route(self):
        """
        Creates a session to find the URL for the loci and schemes
        """
        # Create a new session
        session = OAuth1Session(self.consumer_key,
                                self.consumer_secret,
                                access_token=self.session_token,
                                access_token_secret=self.session_secret)
        # Use the test URL in the GET request
        r = session.get(self.test_rest_url,
                        verify=False)
        if r.status_code == 200 or r.status_code == 201:
            if re.search('json', r.headers['content-type'], flags=0):
                decoded = r.json()
            else:
                decoded = r.text
            # Extract the URLs from the returned data
            self.loci = decoded['loci']
            self.profile = decoded['schemes']

    def download_profile(self):
        """
        Download the profile from the database
        """
        # Set the name of the profile file
        profile_file = os.path.join(self.output_path, 'profile.txt')
        size = 0
        # Ensure that the file exists, and that it is not too small; likely indicating a failed download
        try:
            stats = os.stat(profile_file)
            size = stats.st_size
        except FileNotFoundError:
            pass
        # Only download the profile if the file doesn't exist, or is likely truncated
        if not os.path.isfile(profile_file) or size <= 100:
            # Create a new session
            session = OAuth1Session(self.consumer_key,
                                    self.consumer_secret,
                                    access_token=self.session_token,
                                    access_token_secret=self.session_secret)
            # The profile file is called profiles_csv on the server. Updated the URL appropriately
            r = session.get(self.profile + '/1/profiles_csv',
                            verify=False)
            # On a successful GET request, parse the returned data appropriately
            if r.status_code == 200 or r.status_code == 201:
                if re.search('json', r.headers['content-type'], flags=0):
                    decoded = r.json()
                else:
                    decoded = r.text
                # Write the profile file to disk
                with open(profile_file, 'w') as profile:
                    profile.write(decoded)

    def find_loci(self):
        """
        Finds the URLs for all allele files
        """
        session = OAuth1Session(self.consumer_key,
                                self.consumer_secret,
                                access_token=self.session_token,
                                access_token_secret=self.session_secret)
        # Use the URL for all loci determined above
        r = session.get(self.loci,
                        verify=False)
        if r.status_code == 200 or r.status_code == 201:
            if re.search('json', r.headers['content-type'], flags=0):
                decoded = r.json()
            else:
                decoded = r.text
            # Extract all the URLs in the decoded dictionary under the key 'loci'
            for locus in decoded['loci']:
                # Add each URL to the list
                self.loci_url.append(locus)

    def download_loci(self):
        """
        Uses a multi-threaded approach to download allele files
        """
        # Setup the multiprocessing pool.
        pool = multiprocessing.Pool(processes=self.threads)
        # Map the list of loci URLs to the download method
        pool.map(self.download_threads, self.loci_url)
        pool.close()
        pool.join()

    def download_threads(self, url):
        """
        Download the allele files
        """
        # Set the name of the allele file - split the gene name from the URL
        output_file = os.path.join(self.output_path, '{}.tfa'.format(os.path.split(url)[-1]))
        # Check to see whether the file already exists, and if it is unusually small
        size = 0
        try:
            stats = os.stat(output_file)
            size = stats.st_size
        except FileNotFoundError:
            pass
        # If the file doesn't exist, or is truncated, proceed with the download
        if not os.path.isfile(output_file) or size <= 100:
            # Create a new session
            session = OAuth1Session(self.consumer_key,
                                    self.consumer_secret,
                                    access_token=self.session_token,
                                    access_token_secret=self.session_secret)
            # The allele file on the server is called alleles_fasta. Update the URL appropriately
            r = session.get(url + '/alleles_fasta',
                            verify=False)
            if r.status_code == 200 or r.status_code == 201:
                if re.search('json', r.headers['content-type'], flags=0):
                    decoded = r.json()
                else:
                    decoded = r.text
                # Write the allele to disk
                with open(output_file, 'w') as allele:
                    allele.write(decoded)

    def __init__(self, args):
        self.test_rest_url = 'https://rest.pubmlst.org/db/pubmlst_rmlst_seqdef'
        self.test_web_url = 'https://pubmlst.org/bigsdb/bigsdb.pl?db=pubmlst_rmlst_seqdef'
        self.request_token_url = self.test_rest_url + '/oauth/get_request_token'
        self.session_token_url = self.test_rest_url + '/oauth/get_session_token'
        self.access_token_url = self.test_rest_url + '/oauth/get_access_token'
        self.authorize_url = self.test_web_url + '&page=authorizeClient'
        self.secret_file = args.secret_file
        self.file_path = args.file_path
        self.credentials_file = os.path.join(self.file_path, 'credentials.txt')
        self.credentials_key = os.path.join(self.file_path, 'credentials.key')
        self.output_path = args.output_path
        self.consumer_key = str()
        self.consumer_secret = str()
        self.access_secret = str()
        self.access_token = str()
        self.session_secret = str()
        self.session_token = str()
        self.request_secret = str()
        self.request_token = str()
        self.loci = str()
        self.profile = str()
        self.loci_url = list()
        self.threads = multiprocessing.cpu_count()
