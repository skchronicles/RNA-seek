# Common helper functions shared across the entire workflow

def check_existence(filename):
    """Checks if file exists on filesystem
    :param filename <str>: Name of file to check
    """
    if not os.path.exists(filename):
        sys.exit("File: {} does not exists!".format(filename))


def check_readaccess(filename):
    """Checks permissions to see if user can read a file
    :param filename <str>: Name of file to check
    """
    check_existence(filename)
    if not os.access(filename,os.R_OK):
        sys.exit("File: {} exists, but user cannot read from file due to permissions!".format(filename))


def check_writeaccess(filename):
    """Checks permissions to see if user can write to a file
    :param filename <str>: Name of file to check
    """
    check_existence(filename)
    if not os.access(filename,os.W_OK):
        sys.exit("File: {} exists, but user cannot write to file due to permissions!".format(filename))


def provided(samplelist, condition):
    '''
    Determines if optional rules should run. If an empty list is provided to rule all,
    snakemake will not try to generate that set of target files. If a given condition
    is not met (i.e. False) then it will not try to run that rule.
    '''

    if not condition:
        # If condition is False, returns an empty list to prevent rule from running
        samplelist = []

    return samplelist


def abstract_location(file_address, *args, **kwargs):
    '''
    Determines if a provided file or list of file(s) resides in a remote location.
    If file(s) are determined to reside in remote store, like a S3 or Google Cloud
    Storage, Snakemake's remote wrapper is used to defined remote files.
    This can be extended further to support more file types listed here.
    https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html
    Supported remotes file options include: s3, gs, and sftp
    Input: File path <str> or a list or file paths list[<str>]
    Output: List of files or remote objects
    '''

    # Check if user provided any input
    if not file_address or file_address is None:
        raise IOError("Failed to provide any input files! Input(s) are required to resolve required files.".format(file_address))

    # If given file path to one file, convert it a list[<str>]
    file_list = [file_address] if isinstance(file_address, str) else file_address

    # Loop through list of provided files, and if a remote storage option
    # is given, convert its index to a remote file object.
    for i, uri in enumerate(file_list):
        if uri.lower().startswith('s3://'):
            # Remote option for S3 storage
            import snakemake.remote.S3
            import botocore.session

            if botocore.session.get_session().get_credentials():
                # AWS cli or boto has been configured on target system
                # See ~/.aws/credentials or ~/.boto
                # https://boto.readthedocs.io/en/latest/boto_config_tut.html
                remote_provider = snakemake.remote.S3.RemoteProvider()
            else:
                # If botocore cannot find credentials, try connecting unsigned.
                # This will work for anonymous S3 resources if the resources in the
                # s3 bucket are configured correctly.
                # If a file in provieded as input to a Snakemake rule, only read
                # access is needed to access the remote S3 object.
                remote_provider = snakemake.remote.S3.RemoteProvider(config=botocore.client.Config(signature_version=botocore.UNSIGNED))
            file_list[i] = remote_provider.remote(uri, *args, **kwargs)

        elif uri.lower().startswith('gs://'):
            # Remote option for Google Cloud Storage
            import snakemake.remote.GS
            remote_provider = snakemake.remote.GS.RemoteProvider()
            file_list[i] = remote_provider.remote(uri, *args, **kwargs)

        elif uri.lower().startswith('sftp://'):
            # Remote option for SFTP transfers
            import snakemake.remote.SFTP
            remote_provider = snakemake.remote.SFTP.RemoteProvider()
            file_list[i] = remote_provider.remote(uri, *args, **kwargs)

    return file_list
