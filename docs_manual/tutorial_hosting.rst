.. _tutorial_hosting:

================
Hosting Tutorial
================

This tutorial will describe how to setup SCelVis on your own (web server).
This allows you to host your own scRNA-seq data sets, e.g., as online material for a manuscript (we'd be happy if you cite us) or for sharing data with your colleagues.
We will also shortly demonstrate how to protect your SCelVis server instance with a password.

---------------
Some Background
---------------

Depending on your experience with creating web applications or running web sites, you might ask yourself: "What exactly is a web server?" or "I have a website, how does SCelVis fit in?"
If you're an expert, quickly glance over this section such that you are aware of our nomenclature.

Clients, Servers, and The Web
=============================

.. note::

    You can find out much more about how the web works over at the `Mozilla Website <https://developer.mozilla.org/en-US/docs/Learn/Getting_started_with_the_web/How_the_Web_works>`__.

Overall, the web works as shown in the next figure.
Your browser (the client) sends a HTTP (hypertext transfer protocol) request to the server that answers with an HTTP response.

.. aafig::
    :scale: 80
    :align: center

    +--------+   +--------+
    | Client |   | Server |
    +---+----+   +---+----+
        |            |
        |            |
        |  request   |
        X----------->X
        X            X
        X<-----------X
        |  response  |
        |            |

This is exactly what is happening if you run a web application like SCelVis on your own computer and access it through your browser.
What would be drawbacks when running SCelVis like that on a public web server?

For one, network server programs listen on a server IP (internet protocol) address at a specific port.
The default port for HTTP is 80 and the default for HTTPS (HTTP over an encrypted connection) is 443.
If your URL (uniform resource locator, your server's "web address") begins with ``http://`` and you don't specify a port (e.g., as the port ``8080`` in ``http://example.com:8080``) then your browser will use port ``80`` and if your URL starts with ``https://`` then your browser will use port ``443``.
In most cases, you will want to give your user simple URLs such as ``https://scelvis.yourlab.org`` or ``https://yourlab.org/scelvis`` yet still being able to have more than one scelvis instance running (e.g., one per data set).

Also, SCelVis does not support HTTPS.
Such functionality is usually added via reverse proxies.
What does a reverse proxy do?
To keep the explanation simple: a reverse proxy takes the client's request, forwards it to a *backend server* and then forwards the response of the *backend server* to the client.
On the way forward and backward, the reverse proxy can modify the request and response.
Such modifications may include:

- Wrapping the request and response in an encrypted connection.
  The client uses HTTPS to talk to the reverse proxy while the reverse proxy uses HTTP to talk to the backend server.

- Working as a branching point and forwarding the path ``/dataset-1/`` and everything below to one SCelVis server instance (or particular version) while forwarding ``/dataset-2`` to SCelVis server instance.
  It might still do something else for all other addresses, e.g., display your lab's website.

- Adding password protection to certain parts of your website.

The figure from above can be extended as follows.
Note that one important propery of HTTP is that it is "proxyable" and multiple reverse proxies can be daisy-chained (not shown here).

.. aafig::
    :scale: 80
    :align: center

    +--------+   +-------+-------+   +--------+
    | Client |   | Reverse Proxy |   | Server |
    +---+----+   +-------+-------+   +--------+
        |                |                |
        |                |                |
        |    request     |                |
        X--------------->X    request     |
        X                X--------------->X
        X                X                X
        X                X<---------------X
        X<---------------X    response    |
        |    response    |                |
        |                |                |

Using Docker
============

We will show you how to run a SCelVis web server with the Traefik reverse proxy using Docker.
You can find out more about Docker `on the Docker website <https://docs.docker.com/get-docker/>`__ and more about Traefik `on the traefik website <https://docs.traefik.io/>`__.

Docker is a software that allows you to download **container images** that have software with all dependencies preinstalled.
These images can then be run as **containers** where the container runs programs in an isolated environment.
For example, you might have one Traefik reverse proxy container and two SCelVis containers using the same SCelVis image.
To give your SCelVis servers access to the data to display, you have to explicitely make directories from the **host** (the computer running the containers) visible inside the containers.
Also, network ports have to be explicitely exposed to the host and thus also to the outside network and thus world.
You can also explicitely wire together containers with private networks.
We will describe the Docker commands here in brief and you can find out more details by Googling for them.

We assume that you have Docker already installed and your user can simply execute the  ``docker`` command (otherwise you might have to substitute this with ``sudo docker``).
We assume that you are running this on Linux.

====================
Setup Docker Network
====================

First, we create a new docker network called ``web``.
The containers will communicate on this network but will be isolated from the host otherwise (NB: the ``$`` indicates the start of the prompt, you should enter ``docker network create web`` only on your terminal).

.. code-block:: shell

    $ docker network create web                                                                                                          3bdea3cbe9e52f98f91075cfcd3b01b5f8ae0e15498a30c59c05f85256dcced4

=======================
Setup Traefik Container
=======================

Next, we start the reverse proxy with the Traefik software.
Traefik is configured using command line arguments and via container labels.
We will see more of this below.


.. code-block:: shell

    $ docker run \
        --detach \
        --name traefik \
        --restart unless-stopped \
        --network web \
        --publish 0.0.0.0:80:80 \
        --publish 0.0.0.0:443:443 \
        --volume //var/run/docker.sock:/var/run/docker.sock:ro \
        --volume "$PWD/volumes/traefik/letsencrypt:/letsencrypt" \
        --label traefik.http.middlewares.redirect-to-https.redirectscheme.scheme=https \
        --label 'traefik.http.routers.redirs.rule=hostregexp(`{host:.+}`)' \
        --label traefik.http.routers.redirs.entrypoints=web \
        --label traefik.http.routers.redirs.middlewares=redirect-to-https \
        traefik:v2.0.0-rc3 \
            --log.level=DEBUG \
            --entrypoints.web.address=:80 \
            --entrypoints.websecure.address=:443 \
            --providers.docker \
            --certificatesresolvers.le.acme.email=youremail@yourlab.org \
            --certificatesresolvers.le.acme.tlschallenge=true
        2473377fb83561b183660fadf3b04024bd4c4362aeb37c6307a98a7483b47ee6

This will:

- Start a Docker container with the image ``traefik`` with version ``v2.0.0-rc3`` in the background (``--detach``).
  The image will be downloaded from Docker Hub (https://hub.docker.com).
- The container will be named ``traefik`` and connected to the network ``web``.
- In the case of a crash, the server will be restarted until it is explicitely stopped.
- The ports ``80`` and ``443`` from the container will be mapped to ports ``80`` and ``443`` on the server host (these ports must be free, of course).
  The ``0.0.0.0`` indicates that the ports shall be forwarded for all addresses of the host machine.
- The file (actually a Unix socket) ``/var/run/docker.sock`` from the host will be made available as ``/var/run/docker.sock`` in read-only (``ro``) fashion.
  This is required for Traefik to react if other containers are started or stopped.
- The folder ``volumes/traefik/letsencrypt`` relative to the current directory ``$PWD`` will be made available as ``/letsencrypt`` in the container.
  This is where the SSL certificates from letsencrypt will be stored.
- The ``--label`` arguments make all HTTP requests be forwarded to HTTPS.
- The parameters after ``traefik:v2.0.0-rc3`` are given to the ``traefik`` program itself:
    - Increase log verbosity to ``DEBUG``.
    - Listen on port ``80`` and ``443`` (inside the ``traefik`` container only).
    - React on changes of Docker containers.
    - Obtain an SSL certificate from Letsencrypt.
      **You have to adjust ``youremail@example.com`` above.**
      For Letsencrypt to work, your server must be available under a publically accessible domain name (e.g., ``scelvis.yourlab.org``) and be reachable by the letsencrypt server (in other words: be accessible on the internet).

You can see your container running

.. code-block:: shell

    $ docker ps
    CONTAINER ID        IMAGE                COMMAND                  CREATED             STATUS              PORTS                                         NAMES
    2473377fb835        traefik:v2.0.0-rc3   "/traefik --log.leve…"   40 seconds ago      Up 39 seconds       0.0.0.0:80->80/tcp, 0.0.0.0:443->443/tcp   traefik

And inspect the traefik log output with ``docker logs traefik``

.. code-block:: shell

    $ docker logs traefik
    time="2020-05-26T20:23:04Z" level=info msg="Configuration loaded from flags."
    time="2020-05-26T20:23:04Z" level=info msg="Traefik version 2.0.0-rc3 built on 2019-09-10T17:10:04Z"
    ...


=======================
Setup SCelVis Container
=======================

We will now setup a SCelVis container with the following features:

- Use a public data set in HDF5 format (`download link <https://github.com/bihealth/scelvis/blob/master/examples/hgmm_1k.h5ad?raw=true>`__).
  See :ref:`tutorial_convert` on how to create an appropriate HDF5 file.
- Use a custom Markdown document on the start page that you can modify.
  For example, you can place a link to your lab's website here.
- Allow you to embed custom images in the Markdown document, e.g., your lab's logo, or a picture from your publication.
- Be available below the path `/dataset-1/`.
  This way, you could add a `/dataset-2/` etc.

We will need to perform only a few steps:

1. Prepare a directory with the HDF5 file, your custom Markdown file, and a static image file.
2. Start the SCelVis docker container appropriately.
3. Look at your data.

----------------------------
Preparing the Data Directory
----------------------------

We will create the data directory in your home directory, but any place will do.

.. code-block:: shell

    $ mkdir -p ~/scelvis-data/dataset-1/static
    $ wget -O ~/scelvis-data/dataset-1/hgmm_1k.h5ad \
        'https://github.com/bihealth/scelvis/blob/master/examples/hgmm_1k.h5ad?raw=true'
    $ wget -O ~/scelvis-data/dataset-1/static/rna.png \
        'https://upload.wikimedia.org/wikipedia/commons/a/a4/Pre-mRNA-1ysv-tubes.png'
    $ cat >~/scelvis-data/dataset-1/home.md <<EOF
    ## Your first SCelVis Instance

    Here is an example image:

    ![RNA](/static/rna.png)
    EOF

Look at the resulting structure:

.. code-block:: shell

    $ tree ~/scelvis-data/dataset-1
    /home/user/scelvis-data/dataset-1
    ├── home.md
    ├── hgmm_1k.h5ad
    └── static
        └── rna.png

-----------------------
Start SCelVis Container
-----------------------

Let us fire up a SCelVis Docker container now:


.. code-block:: shell

    $ docker run \
        --detach \
        --name scelvis-dataset-1 \
        --restart unless-stopped \
        --network web \
        --volume $HOME/scelvis-data/dataset-1:/data:ro \
        --env SCELVIS_URL_PREFIX=/dataset-1 \
        --label traefik.enable=true \
        --label 'traefik.http.routers.scelvis-dataset-1.rule=Host(`scelvis.yourlab.org`) && PathPrefix(`/dataset-1`)' \
        --label traefik.http.middlewares.scelvis-dataset-1-stripprefix.stripprefix.prefixes=/dataset-1 \
        --label traefik.http.routers.scelvis-dataset-1.middlewares=scelvis-dataset-1-stripprefix,xforward \
        --label traefik.http.routers.scelvis-dataset-1.entrypoints=websecure \
        --label traefik.http.routers.scelvis-dataset-1.tls.certresolver=le \
        --label traefik.http.middlewares.xforward.headers.customrequestheaders.X-Forwarded-Proto=https \
        --label traefik.http.services.scelvis-dataset-1.loadbalancer.server.port=8050 \
        quay.io/biocontainers/scelvis:0.8.4--py_0 \
        scelvis run \
            --host 0.0.0.0 \
            --port 8050 \
            --data-source /data \
            --custom-home-md /data/home.md \
            --custom-static-folder /data/static \
            --disable-upload \
            --disable-conversion
    a0d0e7395713e529977807e9ee74966f7fe76dfcf8f72d04e2f529e6cbd8ab2e

This will:

- Create a new container named ``scelvis-dataset-1`` using the ``scelvis`` image version ``0.8.4--py_0`` from the ``biocontainers`` repository at the ``quay.io`` server.
- The container will be wired to the network ``web`` and be restarted unless it is explicitely stopped.
- The container is sent to the background (``--detach``).
- Pass ``/home-user/scelvis-data/dataset-1`` from the host to ``/data`` into the container in read-only (``ro``) mode.
- Set the environment varaible ``SCELVIS_URL_PREFIX`` to ``/dataset-1`` inside the container.
- Various labels are attached to the container to communicate with Traefik:

    - Enable traefik.
    - Route the path ``/dataset-1`` on the domain ``scelvis.yourlab.org`` (make sure to change this to your domain name and that ``scelvis.yourlab.org`` actually points to your server).
    - Remove the prefix ``/dataset-1`` when passed to SCelVis by the Traefik reverse proxy.
    - Add the ``X-Forwarded-Proto`` header as HTTPS to SCelVis so it knows that it should create HTTPS URLs.
    - Use HTTPS as the only entry point and use Letsencrypt for creating SSL certificates.
    - Expose the port ``8050`` from within the container to via the ``/dataset-1`` path on the ``scelvis.yourlab.org`` domain.

- Run the SCelVis command line ``scelvis run`` and

    - Run the SCelVis server on all addresses (inside the container) on port ``8050`` (inside the container).
    - Use the data source directory ``/data`` (inside the container) and custom home Markdown file and static folder.
    - Disable the CPU and storage hungry upload and conversion features.

After the container has started, you should be able to navigate to ``https://scelvis.yourlab.org/dataset-1/`` and see your custom SCelVis start page.
In case of problems inspect your container with ``docker logs scelvis-dataset-1`` (and also consider looking into the Traefik logs).

.. code-block::

    $ docker ps
    CONTAINER ID        IMAGE                                       COMMAND                  CREATED             STATUS              PORTS                                      NAMES
    a0d0e7395713        quay.io/biocontainers/scelvis:0.8.4--py_0   "scelvis run --host …"   12 seconds ago      Up 12 seconds                                                  scelvis-dataset-1
    9389b1c23098        traefik:v2.0.0-rc3                          "/traefik --log.leve…"   37 minutes ago      Up 37 minutes       0.0.0.0:80->80/tcp, 0.0.0.0:443->443/tcp   traefik

You can stop the container again:

.. code-block:: shell

    $ docker stop scelvis-dataset-1

Note that you can start as many SCelVis containers as you want (and have memory for).
However, each has to use a different port in the ``web`` network.
Traefik will react on the starting and stopping of the containers by adding new entries into its routes table and thus expose the container to the internet.

==========================
Adding Password Protection
==========================

A quick and easy way to add password protection is to use the Traefik `basicauth middleware <https://docs.traefik.io/middlewares/basicauth/>`__.
For this, you will need to install the ``htpasswd`` program (``sudo apt-get install apache2-utils`` on Ubuntu and ``sudo yum install httpd-tools``).
Use this for creating an htpasswd entry for the user ``user`` and their encrypted password.

.. code-block:: shell

    $ htpasswd -n user
    New password:
    Re-type new password:
    user1:$apr1$sQcbzQ6F$f4frhaAaAOVghxsCnM2Ez/
    $ htpasswd -n user2
    New password:
    Re-type new password:
    user2:$apr1$g.HFrAXQ$HUVZ5mVlOTwWNx/ikTumX1

You can now adjust the ``docker run`` command above (first stop the container with ``docker stop``) to include some more labels to your container:

.. code-block:: shell

    [...] \
    --label "traefik.http.middlewares.scelvis-dataset-1-auth.basicauth.users=user1:$apr1$sQcbzQ6F$f4frhaAaAOVghxsCnM2Ez/,user2:$apr1$g.HFrAXQ$HUVZ5mVlOTwWNx/ikTumX1" \
    --label "traefik.http.middlewares.scelvis-dataset-1-auth.basicauth.realm=Please enter your user name and password to proceed." \
    --label traefik.http.routers.scelvis-dataset-1.middlewares=scelvis-dataset-1-stripprefix,xforward,scelvis-dataset-1-auth \
    [...]

The first two labels are new and define the list of user/(encrypted) password pairs to access the SCelVis server by defining properties of the ``scelvis-dataset-1-auth`` middleware.
The second label is a change to the one above, having the new ``scelvis-dataset-1-auth`` middleware added.
