import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ncbi-analysis",
    version="0.0.1",
    author="Xin Liu",
    author_email="tmeteorj@gmail.com",
    description="An experiment sets of ncbi database",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tmeteorj/ncbi-analysis",
    packages=setuptools.find_packages(),   # 指定需要安装的模块
    install_requires=['peppercorn'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Researchers in Biology',
        'Topic :: Data Mining :: Biology',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        "Operating System :: OS Independent",
    ],
)