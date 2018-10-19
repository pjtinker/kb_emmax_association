FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

# RUN apt-get update
RUN curl -O https://www.cog-genomics.org/static/bin/plink181012/plink_linux_x86_64.zip \
    && unzip plink_linux_x86_64.zip \
    && mv plink /kb/deployment/bin 

RUN wget http://csg.sph.umich.edu//kang/emmax/download/emmax-beta-07Mar2010.tar.gz \
    && tar -xzf emmax-beta-07Mar2010.tar.gz --strip-components=1 \
    && mv -t /kb/deployment/bin emmax emmax-kin 

RUN pip install pandas
# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
