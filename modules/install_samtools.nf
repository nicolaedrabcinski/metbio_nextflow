process SAMTOOLS_INSTALL {
    tag "Installing samtools"
    
    output:
    path "samtools", emit: samtools_bin
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "🧬 Installing samtools..."
    
    # Check if samtools is already available
    if command -v samtools &> /dev/null; then
        echo "✅ samtools found in system PATH"
        mkdir -p samtools
        ln -sf \$(which samtools) samtools/samtools
    else
        # Install samtools via conda if available, otherwise compile from source
        if command -v conda &> /dev/null; then
            echo "📦 Installing samtools via conda..."
            conda create -n samtools_env samtools -c bioconda -y
            
            # Create wrapper script
            mkdir -p samtools
            cat > samtools/samtools << 'EOF'
#!/bin/bash
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate samtools_env
samtools "\$@"
EOF
            chmod +x samtools/samtools
        else
            # Install via package manager or compile from source
            if command -v apt-get &> /dev/null; then
                echo "📦 Installing samtools via apt..."
                sudo apt-get update && sudo apt-get install -y samtools
                mkdir -p samtools
                ln -sf \$(which samtools) samtools/samtools
            elif command -v yum &> /dev/null; then
                echo "📦 Installing samtools via yum..."
                sudo yum install -y samtools
                mkdir -p samtools
                ln -sf \$(which samtools) samtools/samtools
            else
                echo "📦 Compiling samtools from source..."
                # Download and compile samtools
                mkdir -p samtools_build samtools
                cd samtools_build
                
                wget -q https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
                tar -xjf samtools-1.18.tar.bz2
                cd samtools-1.18
                
                # Configure and compile
                ./configure --prefix=\$(pwd)/../../samtools
                make -j\${NPROC:-2}
                make install
                
                cd ../../
                rm -rf samtools_build
            fi
        fi
    fi
    
    # Verify installation
    ./samtools/samtools --version
    echo "✅ samtools installed successfully!"
    """
}