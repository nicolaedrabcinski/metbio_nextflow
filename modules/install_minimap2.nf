process INSTALL_MINIMAP2 {
    tag "Installing minimap2"
    
    output:
    path "minimap2", emit: minimap2_bin
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "🧬 Installing minimap2..."
    
    # Check if minimap2 is already available
    if command -v minimap2 &> /dev/null; then
        echo "✅ minimap2 found in system PATH"
        mkdir -p minimap2
        ln -sf \$(which minimap2) minimap2/minimap2
    else
        # Download and install minimap2
        mkdir -p minimap2
        cd minimap2
        
        # Determine architecture
        ARCH=\$(uname -m)
        if [[ "\$ARCH" == "x86_64" ]]; then
            MINIMAP2_URL="https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2"
            wget -q "\$MINIMAP2_URL" -O minimap2.tar.bz2
            tar -xjf minimap2.tar.bz2 --strip-components=1
            rm minimap2.tar.bz2
        else
            # Compile from source for other architectures
            echo "📦 Compiling minimap2 from source for \$ARCH..."
            wget -q https://github.com/lh3/minimap2/archive/v2.26.tar.gz -O minimap2.tar.gz
            tar -xzf minimap2.tar.gz --strip-components=1
            rm minimap2.tar.gz
            make
        fi
        
        cd ..
    fi
    
    # Verify installation
    ./minimap2/minimap2 --version
    echo "✅ minimap2 installed successfully!"
    """
}